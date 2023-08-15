import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import ttk

# Global parameters (initial values)
L = 10
Nx = 500
Nt = 500
dt = 0.01
V0 = 50.0
x0 = 4.0
sigma = 0.5
hbar = 1.0
mass = 1.0
E = 40.0

# Spatial discretization
dx = L / (Nx - 1)
x = np.linspace(0, L, Nx)

# Create main window
root = tk.Tk()
root.title("Quantum Schrödinger Equation Solver")

# Create figure for plotting
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
ax.set_xlabel("Position")
ax.set_ylabel("Probability Density")
ax.set_title("Time-Dependent Schrödinger Equation Solver")
line, = ax.plot([], [], label="Probability Density")
ax.legend()

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()

# Create toolbar
toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
toolbar.update()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Potential energy function
def potential(x):
    return np.where(np.abs(x - x0) < 0.5, V0, 0)

# Hamiltonian matrix construction
def construct_hamiltonian():
    diag = np.ones(Nx) * (hbar**2 / (2 * mass * dx**2)) + potential(x)
    off_diag = -0.5 * hbar**2 / (mass * dx**2) * np.ones(Nx - 1)
    H = np.diag(diag) + np.diag(off_diag, -1) + np.diag(off_diag, 1)
    return H

# Time evolution using the Crank-Nicolson method
def time_evolution(H, psi0):
    psi = psi0.copy()
    I = np.identity(Nx)
    A = I - 0.5j * dt / hbar * H
    B = I + 0.5j * dt / hbar * H
    for _ in range(Nt):
        psi = np.linalg.solve(A, np.dot(B, psi[:, np.newaxis])).flatten()
    return psi

# Initial wave function (Gaussian wave packet)
def initial_wave_packet(x):
    return np.exp(-0.5 * ((x - L / 2) / sigma)**2) * np.exp(1j * k0 * x)

# Calculate wavenumber k0
potential_energy_x0 = potential(x0)
if E > potential_energy_x0:
    k0 = np.sqrt(2 * mass * (E - potential_energy_x0) / hbar**2)
else:
    k0 = np.sqrt(-2 * mass * (E - potential_energy_x0) / hbar**2) * 1j

# Calculate the time-dependent solution and update the plot
def calculate_and_update_plot():
    global L, Nx, Nt, dt, V0, x0, sigma, hbar, mass, E
    
    # Update parameter values from the Entry widgets
    L = float(param_entries[0].get())
    Nx = int(param_entries[1].get())
    Nt = int(param_entries[2].get())
    dt = float(param_entries[3].get())
    V0 = float(param_entries[4].get())
    x0 = float(param_entries[5].get())
    sigma = float(param_entries[6].get())
    hbar = float(param_entries[7].get())
    mass = float(param_entries[8].get())
    E = float(param_entries[9].get())
    
    H = construct_hamiltonian()
    psi0 = initial_wave_packet(x)
    psi_final = time_evolution(H, psi0)
    prob_density = np.abs(psi_final)**2
    
    ax.clear()
    ax.plot(x, prob_density, label="Probability Density")
    ax.set_xlabel("Position")
    ax.set_ylabel("Probability Density")
    ax.set_title("Time-Dependent Schrödinger Equation Solver")
    ax.legend()
    canvas.draw()

# Create Calculate button
calculate_button = tk.Button(root, text="Calculate and Update Plot", command=calculate_and_update_plot)
calculate_button.pack()

# Create tab control
tab_control = ttk.Notebook(root)
tab_control.pack()

param_tab = ttk.Frame(tab_control)
tab_control.add(param_tab, text="Parameters")

param_frame = tk.Frame(param_tab)
param_frame.pack(padx=20, pady=20)

param_labels = [
    ("L", "Length of the domain"),
    ("Nx", "Number of spatial grid points"),
    ("Nt", "Number of time steps"),
    ("dt", "Time step size"),
    ("V0", "Potential barrier height"),
    ("x0", "Center of the potential barrier"),
    ("sigma", "Width of the wave packet"),
    ("hbar", "Reduced Planck's constant"),
    ("mass", "Particle mass"),
    ("E", "Energy eigenvalue for the bound state")
]
param_entries = []

for i, (label, explanation) in enumerate(param_labels):
    row = tk.Frame(param_frame)
    lab = tk.Label(row, width=15, text=label, anchor='w')
    lab_explanation = tk.Label(row, text=explanation, anchor='w', wraplength=200)
    ent = tk.Entry(row)
    ent.insert(0, str(eval(label)))  # Display the current value
    
    row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
    lab.pack(side=tk.LEFT)
    lab_explanation.pack(side=tk.LEFT)
    ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
    
    param_entries.append(ent)

# Main loop
root.mainloop()
