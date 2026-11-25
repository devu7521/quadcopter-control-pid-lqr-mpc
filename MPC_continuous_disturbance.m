% mpc_quad_v3.m
% MPC quadrotor hover stabilisation with improved terminal cost and constraints.
%
% Main changes versus the previous version:
%   1. Uses DARE terminal cost instead of dlyap on marginally stable Ad.
%   2. Corrects the input-rate penalty gradient.
%   3. Uses a projected fast-gradient QP solver for box and rate constraints.
%   4. Keeps the analytic Cholesky solution as a warm start.
%
% No Optimization Toolbox is required.

clear; close all; clc;

%% ---- Load variables ----------------------------------------
S = load('allvars.mat');
fprintf('Loaded: '); disp(fieldnames(S)');

if isfield(S,'m'),  m = S.m;  else, error('m missing'); end
if isfield(S,'g'),  g = S.g;  else, g = 9.81;           end

if isfield(S,'I1') && isfield(S,'I2') && isfield(S,'I3')
    Ixx = S.I1; Iyy = S.I2; Izz = S.I3;
elseif isfield(S,'I')
    Ixx = S.I;  Iyy = S.I;  Izz = S.I;
else
    error('Need inertia values (I1/I2/I3 or I)');
end

%% ---- Continuous hover linearisation ------------------------
nx = 12; nu = 4;

Ac = zeros(nx);
Ac(1:3, 4:6) = eye(3);
Ac(4, 8)  =  g;          % x-accel couples to pitch
Ac(5, 7)  = -g;          % y-accel couples to roll
Ac(7, 10) = 1;           % phi_dot   = p
Ac(8, 11) = 1;           % theta_dot = q
Ac(9, 12) = 1;           % psi_dot   = r

Bc = zeros(nx, nu);
Bc(6,  1) = 1/m;         % thrust deviation -> z-accel
Bc(10, 2) = 1/Ixx;       % tau_phi
Bc(11, 3) = 1/Iyy;       % tau_theta
Bc(12, 4) = 1/Izz;       % tau_psi

%% ---- Discretisation (ZOH via matrix exponential) ----------
Ts = 0.02;
[Ad, Bd] = zoh_disc(Ac, Bc, Ts, nx, nu);
fprintf('Max |eig(Ad)| = %.6f\n', max(abs(eig(Ad))));

%% ---- Sanity: (Ad,Bd) controllability -----------------------
Co = ctrb(Ad, Bd);
fprintf('Controllability rank = %d / %d\n', rank(Co), nx);

%% ---- MPC horizons ------------------------------------------
Np = 60;
Nc = 20;

%% ---- MPC weights -------------------------------------------
Q = diag([40 40 30 1 1 2 60 60 10 1 1 1]);
R = diag([0.05 0.1 0.1 0.1]);
Rd = diag([0.005 0.01 0.01 0.01]);   % input-rate penalty

%% ---- Terminal weight ---------------------------------------
% Ad has integrators at z=1, so dlyap(Ad',Q) is not appropriate. The
% discrete algebraic Riccati equation gives the infinite-horizon LQR
% terminal cost for the stabilisable pair (Ad,Bd).
try
    [P,~,~] = dare(Ad, Bd, Q, R);
    P = (P + P') / 2;
    terminalCostSource = 'DARE';
catch ME
    warning('DARE failed (%s); using P = 50*Q', ME.message);
    P = 50 * Q;
    terminalCostSource = '50Q fallback';
end
fprintf('Terminal cost source: %s\n', terminalCostSource);

%% ---- Input limits ------------------------------------------
T_hover = m * g;
u_min = [-T_hover*0.9; -8; -8; -8];  % deviations from hover input
u_max = [ T_hover;      8;  8;  8];
du_max = [T_hover*0.3; 4; 4; 4];     % max deviation-input change per step

%% ---- Prediction matrices -----------------------------------
fprintf('Building prediction matrices...\n');
[Phi, Gam] = pred_mats(Ad, Bd, Np, Nc, nx, nu);

%% ---- Precompute QP matrices --------------------------------
fprintf('Precomputing QP data...\n');

Qb = blkdiag(kron(eye(Np-1), Q), P);
Rb = kron(eye(Nc), R);
Rdbar = kron(eye(Nc), Rd);

% D*z = [u_1; u_2-u_1; ...; u_Nc-u_Nc-1].
D = eye(nu*Nc) - diag(ones(nu*(Nc-1),1), -nu);
Rdb = D' * Rdbar * D;

H_qp = 2 * (Gam' * Qb * Gam + Rb + Rdb);
H_qp = (H_qp + H_qp') / 2 + 1e-8 * eye(nu*Nc);

F_precomp = 2 * Gam' * Qb * Phi;
L_chol = chol(H_qp, 'lower');

solverName = 'projected fast gradient';
maxQpIter = 80;
qpTol = 1e-7;
fprintf('QP size: %d variables (%d inputs x %d moves)\n', nu*Nc, nu, Nc);
fprintf('Solver: %s, max iterations: %d\n', solverName, maxQpIter);

%% ---- Simulation --------------------------------------------
Nsim = floor(20 / Ts);

x0 = [0.2; -0.2; 0.1;
      0;    0;    0;
      deg2rad(8); deg2rad(-6); deg2rad(10);
      0.2; -0.1; 0.3];

x_ref = zeros(nx, 1);

Xh = zeros(nx, Nsim+1);  Xh(:,1) = x0;
Uh = zeros(nu, Nsim);
qpIters = zeros(1, Nsim);

u_prev = zeros(nu, 1);
z_warm = zeros(nu*Nc, 1);

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

fprintf('Simulating %d steps (%.1f s)...\n', Nsim, Nsim*Ts);

for k = 1:Nsim
    e = Xh(:,k) - x_ref;

    % Rate term: ||D*z - [u_prev; 0; ...]||_Rdbar^2.
    rate_ref = zeros(nu*Nc, 1);
    rate_ref(1:nu) = u_prev;
    f_qp = F_precomp * e - 2 * D' * Rdbar * rate_ref;

    % Analytic unconstrained solution gives a strong warm start.
    z_uncon = L_chol' \ (L_chol \ (-f_qp));
    z0 = 0.75 * z_warm + 0.25 * z_uncon;

    projectFcn = @(z) project_control_sequence(z, u_min, u_max, u_prev, du_max, nu, Nc);
    [z_opt, qpIters(k)] = solve_projected_qp(H_qp, f_qp, z0, projectFcn, maxQpIter, qpTol);

    u_now = z_opt(1:nu);
    Uh(:,k) = u_now;
    u_prev = u_now;
    z_warm = shift_sequence(z_opt, nu, Nc);

    % Nonlinear plant propagation.
    U_abs = [T_hover; 0; 0; 0] + u_now;
    t0 = (k-1)*Ts;
    t1 = k*Ts;

    [~, Xode] = ode45(@(t,x) quad_nl(t, x, U_abs, m, g, Ixx, Iyy, Izz), ...
                  [t0 t1], Xh(:,k), opts);
    Xh(:,k+1) = Xode(end,:)';
end

fprintf('Done. Mean QP iterations: %.1f\n', mean(qpIters));

%% ---- Plots --------------------------------------------------
Th = (0:Nsim) * Ts;
Uabs = Uh;
Uabs(1,:) = Uabs(1,:) + T_hover;

figure('Name','States','NumberTitle','off','Position',[100 300 900 500]);
labels_pos = {'x (m)','y (m)','z (m)'};
labels_euler = {'\phi (deg)','\theta (deg)','\psi (deg)'};
for i = 1:3
    subplot(2,3,i);
    plot(Th, Xh(i,:), 'LineWidth',1.4); grid on; yline(x_ref(i),'--r');
    xlabel('t (s)'); title(labels_pos{i});
    subplot(2,3,i+3);
    plot(Th, rad2deg(Xh(i+6,:)), 'LineWidth',1.4); grid on; yline(0,'--r');
    xlabel('t (s)'); title(labels_euler{i});
end
sgtitle('Quadrotor States - MPC Hover');

figure('Name','Inputs','NumberTitle','off','Position',[100 100 900 350]);
input_labels = {'Thrust (N)','\tau_\phi (Nm)','\tau_\theta (Nm)','\tau_\psi (Nm)'};
for i = 1:4
    subplot(1,4,i);
    stairs(Th(1:end-1), Uabs(i,:), 'LineWidth',1.2); grid on;
    xlabel('t (s)'); title(input_labels{i});
end
sgtitle('Absolute Control Inputs');

%% ---- Convergence metrics -----------------------------------
pos_err = vecnorm(Xh(1:3,:), 2, 1);
att_err = vecnorm(Xh(7:9,:), 2, 1);
fprintf('\nFinal position error:  %.4f m\n', pos_err(end));
fprintf('Final attitude error:  %.4f rad\n', att_err(end));
fprintf('Settling time (pos<2cm): ');
idx = find(pos_err < 0.02, 1, 'first');
if isempty(idx), fprintf('not reached\n'); else, fprintf('%.2f s\n', Th(idx)); end

save('mpc_sim_results.mat','Th','Xh','Uh','Uabs','Ts','Np','Nc', ...
     'Q','R','Rd','P','terminalCostSource','solverName','qpIters');
fprintf('Saved mpc_sim_results.mat\n');

%% ============================================================
%% LOCAL FUNCTIONS
%% ============================================================

function [Ad, Bd] = zoh_disc(Ac, Bc, Ts, nx, nu)
    M = expm([Ac, Bc; zeros(nu, nx+nu)] * Ts);
    Ad = M(1:nx, 1:nx);
    Bd = M(1:nx, nx+1:end);
end

function [Phi, Gam] = pred_mats(A, B, Np, Nc, nx, nu)
    Phi = zeros(nx*Np, nx);
    Gam = zeros(nx*Np, nu*Nc);

    Markov = cell(Np, 1);
    Atmp = eye(nx);
    for i = 1:Np
        Markov{i} = Atmp * B;
        Atmp = Atmp * A;
    end

    Apow = eye(nx);
    for i = 1:Np
        Apow = A * Apow;
        rows = (i-1)*nx + (1:nx);
        Phi(rows, :) = Apow;

        for j = 1:min(i, Nc)
            cols = (j-1)*nu + (1:nu);
            Gam(rows, cols) = Markov{i-j+1};
        end
    end
end

function [z, iter] = solve_projected_qp(H, f, z0, projectFcn, maxIter, tol)
    z = projectFcn(z0);
    y = z;
    theta = 1;
    L = max(eig(H));
    alpha = 1 / L;

    for iter = 1:maxIter
        grad = H * y + f;
        zNext = projectFcn(y - alpha * grad);

        if norm(zNext - z, inf) < tol
            z = zNext;
            return;
        end

        thetaNext = (1 + sqrt(1 + 4 * theta^2)) / 2;
        y = zNext + ((theta - 1) / thetaNext) * (zNext - z);
        y = projectFcn(y);
        z = zNext;
        theta = thetaNext;
    end
end

function z = project_control_sequence(z, u_min, u_max, u_prev, du_max, nu, Nc)
    z = reshape(z, nu, Nc);
    prev = u_prev;

    for j = 1:Nc
        uj = min(u_max, max(u_min, z(:,j)));
        uj = min(prev + du_max, max(prev - du_max, uj));
        uj = min(u_max, max(u_min, uj));
        z(:,j) = uj;
        prev = uj;
    end

    z = z(:);
end

function zShift = shift_sequence(z, nu, Nc)
    z = reshape(z, nu, Nc);
    zShift = [z(:,2:end), z(:,end)];
    zShift = zShift(:);
end

function xdot = quad_nl(t, x, U, m, g, Ixx, Iyy, Izz)
    phi = x(7); theta = x(8); psi = x(9);
    p = x(10); q = x(11); r = x(12);
    T = U(1); tau = U(2:4);
    if t >= 2 
        Fd = [35; 0; 0];     %  constant X-direction disturbance
    else
        Fd = [0; 0; 0];
    end
    theta = max(-1.4, min(1.4, theta));
    ct = cos(theta);
    if abs(ct) < 1e-3, ct = sign(ct) * 1e-3; end

    cp = cos(phi); sp = sin(phi);
    ct2 = cos(theta); st = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);

    R_mat = [cpsi*ct2, cpsi*st*sp - spsi*cp, cpsi*st*cp + spsi*sp;
             spsi*ct2, spsi*st*sp + cpsi*cp, spsi*st*cp - cpsi*sp;
             -st,      ct2*sp,                ct2*cp];
    acc = (1/m) * (R_mat * [0; 0; T] + Fd) - [0; 0; g];

    E = [1, sp*st/ct, cp*st/ct;
         0, cp,       -sp;
         0, sp/ct,     cp/ct];

    I_vec = [Ixx; Iyy; Izz];
    omega = [p; q; r];
    omega_dot = [1/Ixx; 1/Iyy; 1/Izz] .* (tau - cross(omega, I_vec .* omega));

    xdot = zeros(12,1);
    xdot(1:3) = x(4:6);
    xdot(4:6) = acc;
    xdot(7:9) = E * omega;
    xdot(10:12) = omega_dot;
end