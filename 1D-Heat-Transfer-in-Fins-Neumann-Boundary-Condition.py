import numpy as np
import matplotlib.pyplot as plt

#=======================================================================================
# Input parameters

# Geometrical properties
L = 0.1  # Length of the fin (in meters)
D = 0.005 # Diameter of the fin (in meters)
Nx = 10  # Number of mesh points along the rod
dx = L / (Nx - 1)  # Mesh spacing

# Material and thermal properties
thermal_conductivity = 14.0  # Thermal conductivity for steel (W/m·K)
convection_coef = 100.0 # Convection Coefficient (W/m^2.K)

# Boundary conditions values
T_env = 298.0 # Environment temperature (K)
T_base = 373.0 # Base temperature (K)

#=======================================================================================
# Evaluation of constant value (hP/kA_c)(dx)^2

const = ((convection_coef*(D*np.pi))/(thermal_conductivity*(D/2)**2*np.pi))*(dx)**2

#=======================================================================================
# Define the Gauss-Seidel solver

def gaussSeidel(A, b, x, N, tol): # A: Coefficint (LHS) matrix, B: RHS vector, x: Initial Guess, N: Number of equations, tol: Tolerance
    maxIterations = 1000000
    xprev = [0.0 for i in range(N)]
    for i in range(maxIterations):
        for j in range(N):
            xprev[j] = x[j]
        for j in range(N):
            summ = 0.0
            for k in range(N):
                if (k != j):
                    summ = summ + A[j][k] * x[k]
            x[j] = (b[j] - summ) / A[j][j]
        diff1norm = 0.0
        oldnorm = 0.0
        for j in range(N):
            diff1norm = diff1norm + abs(x[j] - xprev[j])
            oldnorm = oldnorm + abs(xprev[j])  
        if oldnorm == 0.0:
            oldnorm = 1.0
        norm = diff1norm / oldnorm
        if (norm < tol) and i != 0:
            return x

#=======================================================================================
# Solve the problem

n = Nx - 1           # The size of coeeficint matrix

A = np.zeros((n, n))

# Fill the matrix using a for loop and if conditions
for i in range(n):
    for j in range(n):
        if i == j:  # Diagonal elements
            if i == 0:
                A[i, j] = -3 + ((convection_coef / thermal_conductivity)*(2*dx))
            else:
                A[i, j] = -(2 + const)
        elif abs(i - j) == 1:  # Above and below diagonal
            if i == 0 and j == 1:
                A[i, j] = 4    # A_12
            else:
                A[i, j] = 1
        elif i == 0 and j == 2:  # A_13
            A[i, j] = -1
        else:
            A[i, j] = 0
            
B = np.ones(n) * (-const*T_env)  # RHS Vector
B[0] = (2*convection_coef*T_env*dx)/thermal_conductivity
B[-1] = (-const*T_env) - T_base

Guess = np.zeros(n) # Initial Guess

T = gaussSeidel(A, B, Guess, n, 10**-8)
#=======================================================================================
# Show results

# Complete the array for all 10 nodes and reverse the vector
T = T[::-1]
T = np.insert(T, 0, T_base)

# Plot temperature VS location (temperature distribution)
x = np.linspace(0, L, Nx)
plt.plot(x, T)
plt.xlabel("Position along the rod [m]")
plt.ylabel("Temperature [K]")
plt.title("1D Steady-State Heat Transfer Calculated Numerically")
plt.grid(True)
plt.show()

#=======================================================================================
# Compare Numerical and Analytical Solutions

def equation(x): #Calculate real temperatures
    return ((np.cosh(np.sqrt(const/(dx**2))*(L-x)) + (convection_coef/((const/(dx**2))*thermal_conductivity))*np.sinh(np.sqrt(const/(dx**2))*(L-x))) / 
            (np.cosh(np.sqrt(const/(dx**2))*(L)) + (convection_coef/((const/(dx**2))*thermal_conductivity))*np.sinh(np.sqrt(const/(dx**2))*(L))))*(T_base-T_env) + T_env

# Calculate y-axis values using the equation
y_analytical = equation(x)

# Plot the analytical and numerical solutions
plt.plot(x, y_analytical, label='Analytical Solution', color='green', linestyle='--', linewidth=5)
plt.plot(x, T, label='Numerical Solution', color='red')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Analytical and Numerical Solutions")
plt.legend()
plt.grid(True)
plt.show()