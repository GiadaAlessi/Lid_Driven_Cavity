# Numerical Solution of the Lid-Driven Cavity Problem

## Description
This project presents a numerical study of the Lid-Driven Cavity problem, where the steady-state velocity field is computed using the Fractional Step Method. The implemented code calculates the velocity and pressure fields for Reynolds numbers Re = 100 and Re = 400. The solution is obtained using staggered meshes to improve numerical stability, and results are visualized with Gnuplot.

The study includes:
- Implementation of the Fractional Step Method for Navier-Stokes equations.
- Numerical solution of the velocity and pressure fields.
- Analysis of the velocity profiles and flow characteristics.
- Visualization of results using Gnuplot.

## Code Structure
- **`Lid_Driven_Cavity.cpp`** — Main C++ code implementing the convection-diffusion model.
- **`Lid_Driven_Cavity.h`** — Header file containing function declarations and constants.
- **`CMakeLists.txt`** — CMake configuration for building the project.
- **`Es3_Lid_Driven_Cavity_AlessiGiada.pdf`** — Report presenting the results and analysis in detail.

## Results
The detailed results, including velocity profiles for Re = 100 and Re = 400, as well as colormap visualizations, are presented in the report: **`Es3_Lid_Driven_Cavity_AlessiGiada.pdf`**.

## How to Run
1. Ensure **Gnuplot** is installed on your system.
2. Compile the code using CMake:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ./Esercizio4
3. To visualize results, run the Gnuplot script:
   ```bash
   gnuplot VelocityProfilePlot.plt
   gnuplot ColormapPlot.plt

## Author
**Giada Alessi**  
Master in Thermal Engineering  
Universitat Politècnica de Catalunya


