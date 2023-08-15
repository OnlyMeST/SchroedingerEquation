# Schrödinger Equation Solver

The Schrödinger Equation Solver is a Python application that solves the time-dependent Schrödinger equation for a one-dimensional quantum system with a potential barrier. The solver utilizes the Crank-Nicolson method to perform time evolution and visualize the resulting probability density using the matplotlib library.

## Features

- Solves the time-dependent Schrödinger equation for a Gaussian wave packet encountering a potential barrier.
- Uses the Crank-Nicolson method for accurate time evolution.
- Provides a graphical representation of the time-evolving probability density.
- Allows users to customize simulation parameters, such as length, number of grid points, potential barrier properties, and more.

## Prerequisites

Before using the Quantum Schrödinger Equation Solver, ensure you have the following dependencies installed:

- Python 3.x
- NumPy
- Matplotlib
- Tkinter (usually included with Python installations)
- ttk (part of the tkinter module)

You can install the required packages using the following command:

```bash
pip install numpy matplotlib
```

## Usage

1. Clone or download the repository to your local machine.
2. Open a terminal and navigate to the repository's directory.
3. Run the script using the following command:

```bash
python schroedinger.py
```

4. The main application window will appear, showing the initial probability density distribution and customizable parameter fields.
5. Modify the parameters in the "Parameters" tab according to your preferences.
6. Click the "Calculate and Update Plot" button to update the plot based on the new parameters.
7. Observe the time-evolving probability density as the simulation progresses.

## Parameters

The following parameters can be customized using the "Parameters" tab:

- **L**: Length of the domain
- **Nx**: Number of spatial grid points
- **Nt**: Number of time steps
- **dt**: Time step size
- **V0**: Potential barrier height
- **x0**: Center of the potential barrier
- **sigma**: Width of the wave packet
- **hbar**: Reduced Planck's constant
- **mass**: Particle mass
- **E**: Energy eigenvalue for the bound state

Adjust these parameters to explore different scenarios and observe the resulting behavior of the quantum system.

## License

This project is licensed under the MIT License.

## Acknowledgments

The Quantum Schrödinger Equation Solver is based on the Crank-Nicolson method and utilizes the matplotlib and tkinter libraries for visualization and user interface.

## Contact

For questions, suggestions, or issues, please contact me using my [Email](saria.mostafa.pvt@gmail.com).
