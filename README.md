# Quadcopter Flight Control: PID vs LQR vs MPC Under Disturbances

A MATLAB/Simulink simulation of a quadcopter that compares three flight controllers – **PID**, **LQR** and **Model Predictive Control (MPC)** – on how well each one keeps the vehicle stable when different disturbances act on it.

**Author:** Devanshi Yadav · M.Tech, Robotics and Artificial Intelligence, Jaypee Institute of Information Technology (JIIT), Noida

---

## Overview

- **Quadcopter model** built in Simulink (`uav.slx`), with its physical parameters (mass, gravity, moments of inertia) stored in `allvars.mat`.
- **Three controllers** implemented on the same model: PID, LQR and MPC.
- **Disturbance testing:** each controller is run against several disturbance types, and the quadcopter's stability and disturbance rejection are compared.

## Disturbance scenarios

| Scenario | MPC script |
|---|---|
| Constant disturbance | `MPC_constant_disturbance.m` |
| Temporary disturbance | `MPC_temporary_disturbance.m` |
| Continuous disturbance | `MPC_continuous_disturbance.m` |
| Sinusoidal disturbance | `MPC_sinusoidal_disturbance.m` |
| Turbulent wind | `MPC_turbulent_wind.m` |

## Repository contents

| File | Purpose |
|---|---|
| `uav.slx` | Simulink model of the quadcopter |
| `allvars.mat` | Model parameters (mass, gravity, moments of inertia) |
| `LQR_implementation.m` | LQR controller design and simulation, using the parameters in `allvars.mat` |
| `lqr_gains.mat` | Saved LQR gains |
| `lqr_sim_results.mat` | Saved LQR simulation results |
| `MPC_*.m` | MPC simulation under each disturbance scenario (see table above) |
| `mpc_sim_results.mat` | Saved MPC simulation results |

## How to run

1. Open this folder in MATLAB (desktop or MATLAB Online) with Simulink installed. The `lqr` function needs the Control System Toolbox.
2. Run `LQR_implementation.m` to design the LQR controller and simulate it.
3. Run any `MPC_*.m` script to simulate MPC under that disturbance.
4. Open `uav.slx` to explore the Simulink quadcopter model.

## Results

| PID | LQR | MPC (turbulent wind) |
|---|---|---|
| ![PID](results/pid_result.png) | ![LQR](results/lqr_result.png) | ![MPC under turbulent wind](results/mpc_turbulent_wind.png) |

## Skills demonstrated

Quadcopter dynamics and simulation · PID control · Linear Quadratic Regulator (LQR) · Model Predictive Control (MPC) · disturbance rejection and stability analysis · MATLAB · Simulink
