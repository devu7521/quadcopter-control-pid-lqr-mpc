% lqr_from_allvars.m  (fixed: pass parameters into local function)
clear; close all; clc;
% ---- Load your variables ----
S = load('allvars.mat');
vars = fieldnames(S);
disp('Loaded variables:'); disp(vars');

% ---- Map variables (assumptions) ----
if isfield(S,'m');    m = S.m; else error('m (mass) not found in allvars.mat'); end
if isfield(S,'g');    g = S.g; else g = 9.81; warning('g not found, using 9.81'); end
if isfield(S,'I1') && isfield(S,'I2') && isfield(S,'I3')
    Ixx = S.I1; Iyy = S.I2; Izz = S.I3;
elseif isfield(S,'I')
    Ixx = S.I; Iyy = S.I; Izz = S.I;
    warning('Using scalar I for all principal moments of inertia');
else
    error('Inertia variables I1,I2,I3 or I required in allvars.mat');
end

% collect optional K1..K6 if present (not used by default)
Kmotor = [];
for k=1:6
    fname = sprintf('K%d',k);
    if isfield(S,fname)
        Kmotor(end+1) = S.(fname); %#ok<SAGROW>
    end
end
if ~isempty(Kmotor)
    disp('Motor/PID constants (K1..K6) found and stored in Kmotor (not used by default).');
end

% ---- Build linearized state-space (small-angle hover) ----
A = zeros(12);
A(1:3,4:6) = eye(3);
A(4,8) =  g;
A(5,7) = -g;
A(7,10) = 1;
A(8,11) = 1;
A(9,12) = 1;

B = zeros(12,4);
B(6,1)  = 1/m;
B(10,2) = 1/Ixx;
B(11,3) = 1/Iyy;
B(12,4) = 1/Izz;

% ---- LQR design ----
Q = diag([10 10 30 1 1 2 60 60 10 1 1 1]);
R = diag([0.05 0.1 0.1 0.1]);
K = lqr(A,B,Q,R);
disp('Computed LQR gain K (4x12):'); disp(K);
save('lqr_gains.mat','K','m','g','Ixx','Iyy','Izz','Q','R');

% ---- Augmented LQR (integrators on x,y,z) ----
add_integral = true;
if add_integral
    Cpos = [1 0 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0 0 0 0];
    A_aug = [A, zeros(12,3);
             Cpos, zeros(3,3)];
    B_aug = [B; zeros(3,4)];
    Q_aug = blkdiag(Q, diag([100 100 200]));
    R_aug = R;
    K_aug = lqr(A_aug,B_aug,Q_aug,R_aug);
    disp('Computed augmented LQR gain K_aug (4x15).');
    save('lqr_gains.mat','K_aug','-append');
end

% ---- Nonlinear simulation with LQR ----
T_hover = m * g;
x0 = [0.2; -0.2; 0.1; 0; 0; 0; deg2rad(8); deg2rad(-6); deg2rad(10); 0.2; -0.1; 0.3];
tspan = [0 8];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

% Pass the physical parameters into the ODE function explicitly
[t, X] = ode45(@(tt,xx) quad_nl(tt,xx,K,T_hover,m,g,Ixx,Iyy,Izz), tspan, x0, opts);

% compute control history
Uhist = zeros(length(t),4);
for i=1:length(t)
    xnow = X(i,:)';
    delta_u = -K * xnow;
    U = [T_hover;0;0;0] + delta_u;
    Uhist(i,:) = U';
end

% plots
figure('Name','States: positions & Euler','Units','normalized','Position',[0.05 0.05 0.9 0.4]);
subplot(2,3,1); plot(t, X(:,1)); grid on; title('x (m)');
subplot(2,3,2); plot(t, X(:,2)); grid on; title('y (m)');
subplot(2,3,3); plot(t, X(:,3)); grid on; title('z (m)');
subplot(2,3,4); plot(t, rad2deg(X(:,7))); grid on; title('phi (deg)');
subplot(2,3,5); plot(t, rad2deg(X(:,8))); grid on; title('theta (deg)');
subplot(2,3,6); plot(t, rad2deg(X(:,9))); grid on; title('psi (deg)');

figure('Name','Control Inputs','Units','normalized','Position',[0.05 0.05 0.9 0.25]);
subplot(2,2,1); plot(t, Uhist(:,1)); grid on; title('Thrust (N)');
subplot(2,2,2); plot(t, Uhist(:,2)); grid on; title('tau_\phi (N*m)');
subplot(2,2,3); plot(t, Uhist(:,3)); grid on; title('tau_\theta (N*m)');
subplot(2,2,4); plot(t, Uhist(:,4)); grid on; title('tau_\psi (N*m)');

save('lqr_sim_results.mat','t','X','Uhist');

fprintf('Saved lqr_gains.mat and lqr_sim_results.mat.\n');

% ---- Local function: nonlinear dynamics (now receives params explicitly) ----
function xdot = quad_nl(~, x, Klocal, Ttrim, m_local, g_local, Ixx_local, Iyy_local, Izz_local)
    phi = x(7); theta = x(8); psi = x(9);
    p = x(10); q = x(11); r = x(12);
    delta_u = -Klocal * x;
    U = [Ttrim; 0; 0; 0] + delta_u;
    T_total = U(1); tau = U(2:4);
    Rmat = rotmat(phi,theta,psi);
    thrust_body = [0;0;T_total];
    acc = (1/m_local) * (Rmat * thrust_body) - [0;0;g_local];
    I = diag([Ixx_local Iyy_local Izz_local]);
    omega = [p;q;r];
    omega_dot = I \ (tau - cross(omega, I*omega));
    E = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
         0 cos(phi)           -sin(phi);
         0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
    euler_dot = E * omega;
    xdot = zeros(12,1);
    xdot(1:3) = x(4:6);
    xdot(4:6) = acc;
    xdot(7:9) = euler_dot;
    xdot(10:12) = omega_dot;
end

function R = rotmat(phi,theta,psi)
    Rz = [cos(psi) -sin(psi) 0;
          sin(psi)  cos(psi) 0;
             0         0     1];
    Ry = [cos(theta) 0 sin(theta);
             0       1    0;
         -sin(theta) 0 cos(theta)];
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R = Rz * Ry * Rx;
end
