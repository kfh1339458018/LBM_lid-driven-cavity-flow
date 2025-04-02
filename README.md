# LBM_lid-driven-cavity-flow
Lattice Boltzmann Method for simulate lid-driven cavity flow

This is a pratice program when I first learn mesoscopic method. It uses a single relaxation collision model with non-equilibrium extrapolated boundary conditions.
It should be noted that:
1. The basic LBM belongs to an incompressible method and does not consider the energy equation, so the velocity of the top cover should not be set too high.
2. Numerical methods inevitably have numerical errors, and when the grid resolution is low, density is prone to leakage from the four corners. Therefore, it is recommended to have sufficient discrete points (>128) on each edge.

Due to the simplicity of the program, I did not modularize and encapsulate every function. For VSCode users, simply click on CompileRun button to successfully compile and run the C++file.
