#include "Lid_Driven_Cavity.h"
#include<iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;
using namespace std::chrono;

// Function to allocate dynamic matrix
double** AllocateMatrix(int rows, int cols) {
    auto matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols]();
    }
    return matrix;
}

// Function to deallocate dynamic matrix
void DeallocateMatrix(double**& matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    matrix = nullptr;
}

// Function to compute boundary conditions for velocity
void BoundaryConditions(double Uref, int Nx, int Ny, double **u1, double **un, double **un_1, double **v1, double **vn,
                        double **vn_1) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx+1; j++) {
            if (i == 0)
                u1[i][j] = Uref,
                un[i][j] = Uref,
                un_1[i][j] = Uref;
            if (i == Ny - 1 || j == 0 || j == Nx)
                u1[i][j] = 0,
                un[i][j] = 0,
                un_1[i][j] = 0;
        }
    }
    for (int i = 0; i < Ny + 1; i++) {
        for (int j = 0; j < Nx; j++) {
            if (i == 0)
                v1[i][j] = 0,
                vn[i][j] = 0,
                vn_1[i][j] = 0;
            if (i == Ny || j == 0 || j == Nx - 1)
                v1[i][j] = 0,
                vn[i][j] = 0,
                vn_1[i][j] = 0;
        }
    }
}

// Function to compute R(u) on internal nodes of stagg-x mesh
void ComputeRu(double** u, double** v, double** Ru, double rho, double dx, double mu, int Nx, int Ny) {
    for (int i = 1; i < Ny - 1; i++) {
        for (int j = 1; j < Nx; j++) {
            Ru[i][j] = -rho * dx * (1.0 / 2 * (v[i][j - 1] + v[i][j]) * 1.0 / 2 * (u[i][j] + u[i - 1][j])
                                    - 1.0 / 2 * (v[i + 1][j - 1] + v[i + 1][j]) * 1.0 / 2 * (u[i][j] + u[i + 1][j])
                                    + 1.0 / 2 * (u[i][j] + u[i][j + 1]) * 1.0 / 2 * (u[i][j] + u[i][j + 1])
                                    - 1.0 / 2 * (u[i][j] + u[i][j - 1]) * 1.0 / 2 * (u[i][j] + u[i][j - 1]))
                       + mu * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]);
        }
    }
}

//Function to compute R(v) on internal nodes of stagg-y mesh
void ComputeRv (double** u, double** v, double** Rv, double rho, double dx, double mu, int Nx, int Ny) {
    for (int i = 1; i < Ny; i++) {
        for (int j = 1; j < Nx-1; j++) {
            Rv[i][j] = - rho * dx * (1.0 / 2 * (v[i][j] + v[i-1][j]) * 1.0 / 2 * (v[i][j] + v[i-1][j])
                - 1.0 / 2 * (v[i][j] + v[i+1][j]) * 1.0 / 2 * (v[i][j] + v[i+1][j])
                + 1.0 / 2 * (u[i-1][j+1] + u[i][j+1]) * 1.0 / 2 * (v[i][j] + v[i][j+1])
                - 1.0 / 2 * (u[i-1][j] + u[i][j]) * 1.0 / 2 * (v[i][j] + v[i][j-1]))
                + mu * (v[i-1][j] + v[i+1][j] + v[i][j-1] + v[i][j+1] - 4 * v[i][j]);
        }
    }
}

// Function to solve the pressure field with Gauss-Seidel solver on all nodes of stagg-P mesh
void PoissonSolver(double maxResP, double maxIteP, double rho, double dx, double dt, int Nx, int Ny, double **P1,
                   double **vP, double **uP, double **Pg) {
    double resP = maxResP + 1;
    int iteP = 0;
    while (resP > maxResP && iteP < maxIteP) {
        double maxDiffP = 0.0;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (i == 0)
                    P1[i][j] = 1.0 / 1 * (P1[i+1][j] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (j == 0)
                    P1[i][j] = 1.0 / 1 * (P1[i][j+1] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (i == Ny - 1)
                    P1[i][j] = 1.0 / 1 * (P1[i-1][j] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (j == Nx - 1)
                    P1[i][j] = 1.0 / 1 * (P1[i][j-1] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else
                    P1[i][j] = 1.0 / 4 * (P1[i+1][j] + P1[i][j+1] + P1[i-1][j] + P1[i][j-1]
                        - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                double diffP = fabs(P1[i][j] - Pg[i][j]);
                if (diffP > maxDiffP) {
                    maxDiffP = diffP;
                }
            }
        }
        resP = maxDiffP;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                Pg[i][j] = P1[i][j];
            }
        }
        iteP++;
    }
}

// Function to find the interpolated value
double Interpolate(double x, double** data, int dataSize) {
    if (x == data[dataSize - 1][0]) {
        return data[dataSize - 1][1]; // Explicitly handle the last value
    }
    for (int i = 0; i < dataSize - 1; i++) {
        if ((x >= data[i][0] && x <= data[i + 1][0]) ||
            (x <= data[i][0] && x >= data[i + 1][0])) {

            double x1 = data[i][0], x2 = data[i + 1][0];
            double y1 = data[i][1], y2 = data[i + 1][1];
            return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
        }
    }
    return NAN; // If the value is out of range
}

// Function to extract, interpolate the requested values and calculate the mean absolute error
void ResultsManipulation(const string &sourceFile, const string &targetFile, double ref[][3], int refSize,
                               double **data, int N, int index) {
    ifstream inputFile(sourceFile);
    ofstream outputFile(targetFile);

    // Read data from file
    int dataSize = 0;
    double pos, value;
    while (inputFile >> pos >> value && dataSize < N+2) {
        data[dataSize][0] = pos;
        data[dataSize][1] = value;
        dataSize++;
    }

    // Formats output file
    outputFile << left << setw(15) << "Position"
                << setw(20) << "Reference"
                << setw(20) << "Interpolated Result"
                << setw(20) << "Absolute Error"
                << endl;
    outputFile << string(75, '-') << endl;

    // Extraction, interpolation and error calculation
    double sumError = 0.0;
    for (int i = 0; i < refSize; i++) {
        double position = ref[i][0];
        double referenceValue = ref[i][index];
        double interpolatedValue = Interpolate(position, data, dataSize);
        double MAError = fabs(referenceValue - interpolatedValue);

        if (!isnan(interpolatedValue)) {
            outputFile << left << setw(15) << position
                       << setw(20) << referenceValue
                       << setw(20) << interpolatedValue
                       << setw(20) << MAError << endl;
        } else {
            outputFile << left << setw(15) << position
                        << setw(20) << "NAN"
                        << setw(20) << "(out of range)"
                        << endl;
        }
        sumError += MAError;
    }
    outputFile << string(75, '-') << endl;
    outputFile << "Mean Absolute Error = " << sumError / refSize << endl;
    cout << "File " << targetFile << " successfully created." << endl;
    inputFile.close();
    outputFile.close();
}

int main() {
    // Start timer
    auto start = high_resolution_clock::now();

    // Physical data
    int N;
    // Ask mesh refinement
    cout << "Insert number of control volumes: ";
    cin >> N;
    int Nx = N;
    int Ny = N;
    double Re;
    double L = 1.0;
    double dx = L / Nx;
    double rho = 1.0;
    double Uref = 1.0;
    // Ask for desired Reynolds number
    cout << "Insert Reynolds number: ";
    cin >> Re;
    double mu = rho * Uref * L / Re;

    // Numerical data
    // double Cconv = 0.35;
    // double Cdiff = 0.2;
    // CFL condition
    //double dt = min(Cconv * dx / Uref, Cdiff * pow(dx, 2) * rho / mu);    // Works with Re = 400 but not with Re = 100
    double dt = 0.001*0.5;                                                  // Works with both Re = 100 and Re = 400
    double maxRes = 1e-6;
    double maxIte = 1e6;
    double t_count = 0.0;
    double res1 = maxRes + 1;
    double res2 = maxRes + 1;
    int ite = 0;
    double** data = AllocateMatrix(N+2, 2);

    // Pressure, u and v matrices
    double** P1 = AllocateMatrix(Ny, Nx);
    double** Pg = AllocateMatrix(Ny, Nx);
    double** un_1 = AllocateMatrix(Ny, Nx+1);
    double** un = AllocateMatrix(Ny, Nx+1);
    double** u1 = AllocateMatrix(Ny, Nx+1);
    double** uP = AllocateMatrix(Ny, Nx+1);
    double** vn_1 = AllocateMatrix(Ny+1, Nx);
    double** vn = AllocateMatrix(Ny+1, Nx);
    double** v1 = AllocateMatrix(Ny+1, Nx);
    double** vP = AllocateMatrix(Ny+1, Nx);

    // Convective-diffusive term R() matrices
    double** Ru1 = AllocateMatrix(Ny, Nx+1);
    double** Run = AllocateMatrix(Ny, Nx+1);
    double** Run_1 = AllocateMatrix(Ny, Nx+1);
    double** Rv1 = AllocateMatrix(Ny+1, Nx);
    double** Rvn = AllocateMatrix(Ny+1, Nx);
    double** Rvn_1 = AllocateMatrix(Ny+1, Nx);

    // Initial velocity fields = 0 for internal nodes
    for (int i = 1; i < Ny - 1; i++) {
        for (int j = 1; j < Nx; j++) {
            un[i][j] = 0;
            un_1[i][j] = un[i][j];
        }
    }
    for (int i = 1; i < Ny; i++) {
        for (int j = 1; j < Nx - 1; j++) {
            vn[i][j] = 0;
            vn_1[i][j] = vn[i][j];
        }
    }

    // Initial pressure field
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            Pg[i][j] = 0;
        }
    }

    // Apply constant boundary conditions of velocity
    BoundaryConditions(Uref, Nx, Ny, u1, un, un_1, v1, vn, vn_1);

    // Compute R()^n-1
    ComputeRu(un_1, vn_1, Run_1, rho, dx, mu, Nx, Ny);
    ComputeRv(un_1, vn_1, Rvn_1, rho, dx, mu, Nx, Ny);

    // Compute R()^n
    ComputeRu(un, vn, Run, rho, dx, mu, Nx, Ny);
    ComputeRv(un, vn, Rvn, rho, dx, mu, Nx, Ny);

    // Time loop until convergence (aka steady state) is reached
    while (res1 > maxRes && res2 > maxRes && ite < maxIte) {
        double maxDiff1 = 0.0;
        double maxDiff2 = 0.0;

        // Step 1a
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                uP[i][j] = un[i][j] + dt / (rho * pow(dx, 2)) * (3.0 / 2 * Run[i][j] - 1.0 / 2 * Run_1[i][j]);
            }
        }

        // Step 1b
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                vP[i][j] = vn[i][j] + dt / (rho * pow(dx, 2)) * (3.0 / 2 * Rvn[i][j] - 1.0 / 2 * Rvn_1[i][j]);
            }
        }

        // Step 2
        PoissonSolver(maxRes, maxIte, rho, dx, dt, Nx, Ny, P1, vP, uP, Pg);

        // Step 3
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                u1[i][j] = uP[i][j] - dt / rho * (P1[i][j] - P1[i][j - 1]) / dx;
                double diff1 = fabs(u1[i][j] - un[i][j]);
                if (diff1 > maxDiff1) {
                    maxDiff1 = diff1;
                }
            }
        }
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                v1[i][j] = vP[i][j] - dt / rho * (P1[i - 1][j] - P1[i][j]) / dx;
                double diff2 = fabs(v1[i][j] - vn[i][j]);
                if (diff2 > maxDiff2) {
                    maxDiff2 = diff2;
                }
            }
        }
        t_count += dt;
        res1 = maxDiff1;
        res2 = maxDiff2;

        // Update u^n and u^n-1
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                un_1[i][j] = un[i][j];
                un[i][j] = u1[i][j];
            }
        }
        // Update v^n and v^n-1
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                vn_1[i][j] = vn[i][j];
                vn[i][j] = v1[i][j];
            }
        }

        // Compute R()^n+1
        ComputeRu(u1, v1, Ru1, rho, dx, mu, Nx, Ny);
        ComputeRv(u1, v1, Rv1, rho, dx, mu, Nx, Ny);

        // Update R()^n and R()^n-1
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                Run_1[i][j] = Run[i][j];
                Run[i][j] = Ru1[i][j];
            }
        }
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rvn_1[i][j] = Rvn[i][j];
                Rvn[i][j] = Rv1[i][j];
            }
        }

        // Go to next time step, if needed
        ite++;
    }

    cout << "Steady state reached in " << t_count << " seconds" << endl;

    // File to check the meshes
    ofstream TestFile("Test.txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            TestFile << u1[i][j] << " ";
        }
        TestFile << endl;
    }
    TestFile.close();

    // File to print u(y) at L/2 for Plot 1
    ofstream Uyl2("Uy.txt");
    Uyl2 << 1 << " " << 1 << endl;
    for (int i = 0; i < Ny; i++) {
        Uyl2 << L - i * dx - dx / 2 << " " << (u1[i][Nx / 2] + u1[i][Nx / 2 + 1]) / 2 << endl;
    }
    Uyl2 << 0 << " " << 0 << endl;
    Uyl2.close();

    // File to print data for Gnuplot 1
    ofstream GnuplotData1("GnuplotData1.txt");
    double ref1[17][3] = {
        {1.0000, 1.00000, 1.00000},
        {0.9766, 0.84123, 0.75837},
        {0.9688, 0.78871, 0.68439},
        {0.9609, 0.73722, 0.61756},
        {0.9531, 0.68717, 0.55892},
        {0.8516, 0.23151, 0.29093},
        {0.7344, 0.00332, 0.16256},
        {0.6172, -0.13641, 0.02135},
        {0.5000, -0.20581, -0.11477},
        {0.4531, -0.21090, -0.17119},
        {0.2813, -0.15662, -0.32726},
        {0.1719, -0.10150, -0.24299},
        {0.1016, -0.06434, -0.14612},
        {0.0703, -0.04775, -0.10338},
        {0.0625, -0.04192, -0.09266},
        {0.0547, -0.03717, -0.08186},
        {0.0000, 0.00000, 0.00000}
    };

    int index1 = 0;
    if (Re == 100)
        index1 = 1;
    else if (Re == 400)
        index1 = 2;

    for (auto & i : ref1) {
        GnuplotData1 << i[0] << " " << i[index1] << std::endl;
    }
    GnuplotData1.close();

    // Generate Gnuplot script 1
    ofstream Gnuplot1("ComparisonUy.plt");

    Gnuplot1 << "set terminal pngcairo size 800,600 enhanced font 'Arial,14'\n";
    Gnuplot1 << "set output 'Comparison_Uy.png'\n";
    Gnuplot1 << "set xlabel 'u'\n";
    Gnuplot1 << "set ylabel 'y'\n";
    Gnuplot1 << "set title 'Velocity Profile u(y) at L/2 with Re = " << Re << "'\n";
    Gnuplot1 << "set key top left\n";
    Gnuplot1 << "set grid\n";
    Gnuplot1 << "set xrange [-0.4:1]\n";
    Gnuplot1 << "set yrange [0:1]\n";
    Gnuplot1 << "plot 'Uy.txt' using 2:1 with lines lw 2 lc rgb 'blue' title 'Result', "
            << "'GnuplotData1.txt' using 2:1 with points pt 7 ps 1.2 lc rgb 'red' title 'Reference'\n";
    Gnuplot1.close();

    // Automatically run Gnuplot 1
    system("gnuplot ComparisonUy.plt");
    cout << "Gnuplot script 1 executed: 'Comparison_Uy.png' generated." << endl;

    // File to print v(x) at L/2 for Plot 2
    ofstream Vxl2("Vx.txt");
    Vxl2 << 0 << " " << 0 << endl;
    for (int j = 0; j < Nx; j++) {
        Vxl2 << j * dx + dx / 2 << " " << (v1[Ny / 2][j] + v1[Ny / 2 + 1][j]) / 2 << endl;
    }
    Vxl2 << 1 << " " << 0 << endl;
    Vxl2.close();

    // File to print data for Gnuplot 2
    ofstream GnuplotData2("GnuplotData2.txt");
    double ref2[17][3] = {
        {0.0000, 0.00000, 0.00000},
        {0.0625, 0.09233, 0.18360},
        {0.0703, 0.10091, 0.19713},
        {0.0781, 0.10890, 0.20920},
        {0.0938, 0.12317, 0.22965},
        {0.1563, 0.16077, 0.28124},
        {0.2266, 0.17507, 0.30203},
        {0.2344, 0.17527, 0.30174},
        {0.5000, 0.05454, 0.05186},
        {0.8047, -0.24533, -0.38598},
        {0.8594, -0.22445, -0.44993},
        {0.9063, -0.16914, -0.23827},
        {0.9453, -0.10313, -0.22847},
        {0.9531, -0.08864, -0.19254},
        {0.9609, -0.07391, -0.15663},
        {0.9688, -0.05906, -0.12146},
        {1.0000, 0.00000, 0.00000}
    };

    int index2 = 0;
    if (Re == 100)
        index2 = 1;
    else if (Re == 400)
        index2 = 2;

    for (auto & i : ref2) {
        GnuplotData2 << i[0] << " " << i[index2] << std::endl;
    }
    GnuplotData2.close();

    // Generate Gnuplot script 2
    ofstream Gnuplot2("ComparisonVx.plt");

    Gnuplot2 << "set terminal pngcairo size 800,600 enhanced font 'Arial,14'\n";
    Gnuplot2 << "set output 'Comparison_Vx.png'\n";
    Gnuplot2 << "set xlabel 'x'\n";
    Gnuplot2 << "set ylabel 'v(x)'\n";
    Gnuplot2 << "set title 'Velocity Profile v(x) at L/2 with Re = " << Re << "'\n";
    Gnuplot2 << "set key top left\n";
    Gnuplot2 << "set grid\n";
    Gnuplot2 << "set xrange [0:1]\n";
    Gnuplot2 << "set yrange [-0.5:0.4]\n";
    Gnuplot2 << "plot 'Vx.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Result', "
            << "'GnuplotData2.txt' using 1:2 with points pt 7 ps 1.2 lc rgb 'red' title 'Reference'\n";
    Gnuplot2.close();

    // Automatically run Gnuplot 2
    system("gnuplot ComparisonVx.plt");
    cout << "Gnuplot script 2 executed: 'Comparison_Vx.png' generated." << endl;

    // File to print velocity distribution for Plot 3
    ofstream VelocityDistribution("VelocityDistribution.txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            VelocityDistribution << j * dx + dx / 2 << " " << L - i * dx - dx / 2 << " " << sqrt(
                pow(v1[i][j], 2) + pow(u1[i][j], 2)) << endl;
        }
        VelocityDistribution << "\n";
    }
    VelocityDistribution.close();

    // Generate Gnuplot script 3
    ofstream Gnuplot3("MagnitudeMap.plt");

    Gnuplot3 << "set terminal pngcairo size 600,600 enhanced font 'Arial,12'\n";
    Gnuplot3 << "set output 'MagnitudeMap_Plot.png'\n";
    Gnuplot3 << "set pm3d map\n";
    Gnuplot3 << "set palette defined ("
             << "0 'blue', "
             << "0.3 'cyan', "
             << "0.5 'green', "
             << "0.7 'yellow', "
             << "1 'red')\n";
    Gnuplot3 << "set colorbox vertical\n";
    Gnuplot3 << "set xlabel 'X (m)'\n";
    Gnuplot3 << "set ylabel 'Y (m)'\n";
    Gnuplot3 << "set title 'Velocity Magnitude Map with Re = " << Re << "' font 'Arial,20'\n";
    Gnuplot3 << "set xrange [0:1]\n";
    Gnuplot3 << "set yrange [0:1]\n";
    Gnuplot3 << "set autoscale\n";
    Gnuplot3 << "set cbrange [0:1]\n";
    Gnuplot3 << "set size square\n";
    Gnuplot3 << "splot 'VelocityDistribution.txt' using 1:2:3 with pm3d notitle\n";

    Gnuplot3.close();

    // Automatically run Gnuplot 3
    system("gnuplot MagnitudeMap.plt");
    cout << "Gnuplot script 3 executed: 'MagnitudeMap_Plot.png' generated." << endl;

    // Extract and interpolate data for error analysis
    ResultsManipulation("Uy.txt", "FinalResults_Uy.txt", ref1, 17, data, N, index1);
    ResultsManipulation("Vx.txt", "FinalResults_Vx.txt", ref2, 17, data, N, index2);

    // Stop timer and print total duration
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Total Execution Time: " << duration.count() << " seconds" << endl;

    // Deallocate matrix
    DeallocateMatrix(P1, Ny);
    DeallocateMatrix(Pg, Ny);
    DeallocateMatrix(un_1, Ny);
    DeallocateMatrix(un, Ny);
    DeallocateMatrix(u1, Ny);
    DeallocateMatrix(uP, Ny);
    DeallocateMatrix(vn_1, Ny+1);
    DeallocateMatrix(vn, Ny+1);
    DeallocateMatrix(v1, Ny+1);
    DeallocateMatrix(vP, Ny+1);
    DeallocateMatrix(Ru1, Ny);
    DeallocateMatrix(Run, Ny);
    DeallocateMatrix(Run_1, Ny);
    DeallocateMatrix(Rv1, Ny+1);
    DeallocateMatrix(Rvn, Ny+1);
    DeallocateMatrix(Rvn_1, Ny+1);
    DeallocateMatrix(data, N+2);

    return 0;
}