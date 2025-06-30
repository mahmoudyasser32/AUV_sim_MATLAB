# AUV_sim_MATLAB

This repository contains a MATLAB-based simulation environment for Autonomous Underwater Vehicle (AUV) control and tracking. The simulation is built using the open-source [SamSim AUV Simulator](https://github.com/SamMans/SamSim) developed by **Dr. Samer A. Mohamed**, and extended with custom control algorithms developed from scratch.

## 📌 Overview

This project simulates the 6-DOF dynamics of an AUV navigating in a 3D underwater environment. The main contribution of this repository is the implementation of a hybrid **Sliding Mode Control (SMC)** and **Model Predictive Control (MPC)** strategy, designed and tuned specifically for this simulation.

While the simulation environment, plant model, and basic utilities are adapted from SamSim, the control logic is fully original and developed independently.

## 🎯 Features

- Full 6-DOF nonlinear AUV dynamics simulation.
- Hybrid **SMC-MPC controller** designed from scratch.
- Trajectory tracking and position regulation tests.
- Modular structure to allow future controller development or integration.
- Visualization of AUV behavior using MATLAB plotting tools.

## 🛠️ Structure

AUV_sim_MATLAB/
│
├── main.m # Main simulation script
├── AUV_parameters.m # AUV model and parameters
├── dynamics.m # AUV dynamic model
├── controller/
│ ├── Control_SMC_MPC.m # Custom controller implementation (SMC + MPC)
│ └── utils/ # Helper functions for control
├── path_planning/ # Trajectory and reference generation
├── results/ # Stored results and plots
└── README.md # This file

Copy
Edit

## 🔧 Dependencies

- MATLAB (tested with R2021a and newer)
- Optimization Toolbox (for MPC)
- Control System Toolbox

## 📈 Control Approach

The **Control_SMC_MPC** controller combines the robustness of **Sliding Mode Control (SMC)** with the prediction and constraint-handling abilities of **Model Predictive Control (MPC)**. This hybrid control structure ensures smooth trajectory tracking even under disturbances and model uncertainties. The control structure works in two stages:

1. **SMC** handles the fast nonlinear dynamics and introduces robustness.
2. **MPC** plans optimal control actions over a prediction horizon based on linearized or approximated dynamics.

Both components are implemented from scratch without using high-level MPC toolboxes.

## 📚 Credits

- **Simulation base:** [SamSim AUV Simulator](https://github.com/SamMans) by **Dr. Samer A. Mohamed**  
- **Control system design and implementation:** [Mahmoud Yasser](https://github.com/mahmoudyasser32)

## 📬 Contact

If you have any questions or suggestions, feel free to open an issue or contact me via [GitHub](https://github.com/mahmoudyasser32).

---