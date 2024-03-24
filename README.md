# ♨ 1D Heat Transfer in Fins ♨
🟢 Python script to solve the 1D heat equation and gain temperature distribution in a fin with Dirichlet or Neumann boundary condition at tip, using TDMA algorithm.  
🟢 This solution is based on finite difference method.  
🟢 Compare numerical solution with analytical solution.  

  
## 🧬 Input variables to define the problem:  
L: Length of the fin (in meters)  
D: diameter of the fin (in meters)  
Nx: Number of nodes along the rod  

thermal_conductivity: Thermal conductivity for steel (in W/m·K)  
convection_coef: Convection Coefficient (in W/m^2.K)  

T_env: Environment temperature (in K)  
T_base: Base temperature (in K)  
T_tip: Tip temperature (in K) (for Dirichlet B.C only)  

    
## 🤖 Usage  
Just replace the existing values of input parameters with your desired values and run the code.
