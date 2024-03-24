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
T_tip = 308.0 # Tip temperature (K)

#=======================================================================================
# Evaluation of constant value (hP/kA_c)(dx)^2

const = ((convection_coef*(D*np.pi))/(thermal_conductivity*(D/2)**2*np.pi))*(dx)**2

#=======================================================================================
# Define the TDMA solver

def TDMA(a,b,c,d):
    n = len(a)
    ϕ = np.empty(n)
    P = np.empty(n)
    Q = np.empty(n)
    
    P[0] = -c[0]/b[0]
    Q[0] = d[0]/b[0]
    
    #Forward elimination
    for i in range (1, n):
        P[i] = - c[i] / (b[i] + a[i]*P[i-1])
        Q[i] = (d[i] - a[i]*Q[i-1])/(b[i] + a[i]*P[i-1])
    
    # Last point temperature
    ϕ[-1] = Q[-1]
    
    # Backward substitution
    for i in range(n-2, -1, -1):
        ϕ[i] = P[i]*ϕ[i+1] + Q[i]

    return ϕ

#=======================================================================================
# Solve the problem

n = Nx - 2                           # The size of coeeficint matrix

b = np.ones(n) * (-(2+const))        # Diagonal Elements of LHS

a = np.ones(n)                       # Above diagonal elements of LHS
a[0] = 0

c = np.ones(n)                       # Below diagonal elements of LHS
c[-1] = 0 

d = np.ones(n) * (-const*T_env)      # RHS Vector
d[0] = (-const*T_env) - T_base
d[-1] = (-const*T_env) - T_tip

T = TDMA(a,b,c,d)

#=======================================================================================
# Show results

# Complete the array for all 10 nodes
T = np.insert(T, 0, T_base)
T = np.append(T, T_tip)

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
    return ((((T_tip-T_env)/(T_base-T_env))*np.sinh(np.sqrt(const/(dx**2))*x) + np.sinh(np.sqrt(const/(dx**2))*(L-x)))/(np.sinh(np.sqrt(const/(dx**2))*(L))))*(T_base-T_env) + T_env
            
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