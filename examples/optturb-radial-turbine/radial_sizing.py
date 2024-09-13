import numpy as np
import matplotlib.pyplot as plt
 
Omega = 60000             # [rev/min]
m_dot = 1                 # [kg/s]  
P_tin = 1.5e6             # [Pa]
T_tin = 1100              # [K]
P_3 = 0.8e6               # [Pa]
Power = 70000             # [W]
D_max = 3.9               # [in]
L_max = 3.25              # [in]
Clearance_Ratio = 0.02    # Tip clearance at exit (Station 3) (1-H_3/Span_3)
work = Power/m_dot        # [m^2/s^2]
 
# Fluid properties of HeXe (Molecular Weight assumed 40g/mol)
gamma = 1.667
MW = 40
R = 8.314/(MW/1000)
Cp = 20.785/(MW/1000)
 
# Using Whitield"s (1990) design method for sizing the turbine rotor
S = work/(Cp*T_tin)
Z = np.linspace(5,25)
beta_2 = np.arccos(1-2/Z)                               # Relative inlet flow incidence (IFR turbine blade is radial) [Radians]
alpha_2 = 90 - beta_2*180/(np.pi*2)                     # [Degrees]
a_01 = np.sqrt(gamma*R*T_tin)
U_2 = a_01*np.sqrt(S/((gamma-1)*np.cos(beta_2)))        # Rotor tip speed [m/s]
c_theta2 = U_2*(1-2/Z)                                  # [m/s]
c_m2 = c_theta2/np.tan(alpha_2*np.pi/180)               # [m/s]
D_2 = 60*U_2/(np.pi*Omega)                              # Inlet diameter [m]
r_2 = D_2/2
#b_2 = m_dot/(4*np.pi*rho_2*c_m2*r_2**2)                 # Inlet vane height [m]
c_o = U_2/0.707                                         # Spouting velocity [m/s], Based on Rodgers and Geiser (1987) for maximum efficiency
 
# Exit parameters (Station 3)
r_3s = 0.7*r_2                                          # Based on Rohlik's recommendation (1968)
r_3h = 0.4*r_3s                                         # Based on Rohlik's recommendation (1968)
t_3 = Clearance_Ratio*(r_3s-r_3h)
r = np.linspace(r_3h,r_3s-t_3)
c_m3 = 0.25*U_2                                         # Exit meridional velocity [m/s], Based on Rodgers and Geiser (1987) for maximum efficiency
c_3 = c_m3                                              # Assume c_theta3 = 0 (flow is entirely axial at exit)
beta_3 = np.arctan(U_2*r/(c_m3*r_2))*180/np.pi          # Exit blade angle, as a function of radius
 
# Specific Speed Calculation
A_d = np.pi*D_2**2/4                                    # Rotor disk area [m^2]
A_3 = A_d*0.41                                          # Exit area [m^2], Based on Rohlik (1968)
Omega_s = (2*np.sqrt(2))**1.5*(U_2/c_o)**1.5*np.sqrt(c_3*A_3/(Omega*np.pi/60*D_2**3))   # Specific speed
 
eta_ts = S/(1-(P_3/P_tin)**((gamma-1)/gamma))
eta_tt = 1/(1/eta_ts-(c_3**2/(2*work)))
 
print("Inlet Diameter[in] = ",D_2*39.3701)
 
# Efficiency as function of number of blades
plt.figure(1)
plt.plot(Z,Omega_s,color='black')
plt.grid(True)
plt.xlabel("Number of Blades")
plt.ylabel("Total-to-Total Efficiency")
 
# Exit blade angle for each number of blades
plt.figure(2)
plt.plot(r,beta_3)
plt.grid(True)
plt.xlabel("radius [m]")
plt.ylabel("Blade Exit Angle[deg]")
 
# Exit blade angle for each number of blades
plt.figure(3)
plt.plot(Z,D_2*39.3701)
plt.grid(True)
plt.xlabel("Number of Blades")
plt.ylabel("Inlet Radius[in]")
 
# Exit blade angle for each number of blades
plt.figure(4)
plt.plot(Z,Z*L_max/(D_2*39.3701))
plt.grid(True)
plt.xlabel("Number of Blades")
plt.ylabel("Maximum Solidity Parameter")
plt.show()