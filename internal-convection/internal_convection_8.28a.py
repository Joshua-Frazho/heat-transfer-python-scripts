# INTERNAL CONVECTION: Bergman et al. Problem 8.28a
# Joshua Frazho
# 2026/05/04

import numpy as np

# FUNCTIONS
# T_avg:    returns the average of two temperatures
def T_avg(T1, T2):

    return np.average([T1, T2])

# T_bounds: returns the two neighboring numbers of a certain step given a temperature
#           the step value is determined by the temperature interval printed in the corresponding properties table
def T_bounds(T):
    p_step = 50 #K, temperature interval in properties table for air

    T_lower = np.floor(T / p_step) * p_step #K, lower neighboring temperature value found in properties table
    T_upper = np.ceil(T / p_step) * p_step  #K, upper neighboring temperature value found in properties table

    return T_lower, T_upper

# nu_avg_turb:  Dittus-Boelter equation
#               internal convection correlation for the local Nusselt number in a circular tube characterized by fully-developed, turbulent flow
#               for fully-developed flow, the local Nusselt number is assumed to be constant, hence, the correlation for the local Nusselt number can be used for the average Nusselt number
def nu_avg_turb(Re, Pr, n):
    # Nu_avg_D = 0.023*Re_D^(4/5)*Pr^n
    # n = 0.4 for Ts > Tm, n = 0.3 for Ts < Tm
    # Conditions: 0.6 <≈ Pr <≈ 160, Re_D >≈ 1e4, L/D >≈ 10
    nu_avg_t = 0.023 * (Re ** (4/5)) * (Pr ** n)

    return nu_avg_t

# L_required:   returns the required circular tube length to meet each condition from rearranging the outlet temperature equation 
def L_required(Ts, Tmi, Tmo, h, P, m_dot, cp):
    # Rearranged outlet temp equation, solving for L
    # L = -m_dot*cp / (h*P) * ln((Ts - Tmo)/(Ts - Tmi))
    L = -(m_dot * cp) / (h * P) * np.log((Ts - Tmo) / (Ts - Tmi))

    return L

# KNOWN
# Air flows through copper tubes submerged in an ice/water bath
D     = 0.05         #m, tube inner diameter
Ts    = 0 + 273.15   #K, surface temperature (ice/water bath at 0°C)
T_mi  = 24 + 273.15  #K, air inlet temperature (24°C)
T_mo  = 14 + 273.15  #K, air outlet temperature (14°C)
m_dot = 0.01         #kg/s, mass flow rate per tube

A_c = np.pi / 4 * D**2  #m^2, tube cross-sectional area
P   = np.pi * D         #m, wetted perimeter (circular tube)

# Ice tank properties
N_tubes  = 10       #number of tubes
V_tank   = 10       #m^3, total tank volume
ice_frac = 0.80     #initial volume fraction of ice
rho_ice  = 920      #kg/m^3, density of ice
h_sf     = 3.34e5   #J/kg, latent heat of fusion of ice

# OBJECTIVE
# Find tube length L for T_mo = 14°C
# Find time t to completely melt the ice

# SOLUTION
# ASSUMPTIONS
# - steady-state
# - uniform surface temperature (Ts = 0°C, isothermal ice bath)
# - fully-developed flow
# - air properties evaluated at mean temperature Tm = (Tmi + Tmo)/2
# - 10 identical tubes operating in parallel

# ANALYSIS
print("INTERNAL CONVECTION: Problem 8.28a")

# Mean air temperature for property evaluation
T_m = T_avg(T_mi, T_mo)
print(f"Mean air temperature: {T_m:.2f} K")

# Air properties at Tm = 292 K (19°C), from Table A.4
# Interpolating between 250 K and 300 K
T_lower, T_upper = T_bounds(T_m)
print(f"Bounds for Property Table: {T_lower}, {T_upper} K")

# Table A.4: Thermophysical Properties of Gases at Atmospheric Pressure (Air)
# Table A.4 from Bergman, T. L., and Lavine, A. S., 2017, Fundamentals of Heat and Mass Transfer, John Wiley and Sons, Inc.
# Note: it would be ideal to have a Pandas dataframe that can be referenced across multiple scripts so thermophysical properties would not have to be hard-coded each time
rho = np.interp(T_m, [T_lower, T_upper], [1.3947, 1.1614])      #kg/m^3, density
cp  = np.interp(T_m, [T_lower, T_upper], [1006, 1007])          #J/kg·K, specific heat
mu  = np.interp(T_m, [T_lower, T_upper], [159.6e-7, 184.6e-7])  #Ns/m^2, dynamic viscosity
k   = np.interp(T_m, [T_lower, T_upper], [22.3e-3, 26.3e-3])    #W/mK, thermal conductivity
Pr  = np.interp(T_m, [T_lower, T_upper], [0.72, 0.707])         #unitless, Prandtl number 

# Calculating the Reynolds number to check the flow regime
u_m  = m_dot / (rho * A_c)  #m/s, mean velocity
Re_D = (u_m * rho * D) / mu #Reynolds number

print(f"Mean velocity: {u_m:.3f} m/s")
print(f"Reynolds number: {Re_D:.4e}")

if Re_D < 2300:
    print(f"Laminar flow: {Re_D:.2e} < 2300")
elif Re_D > 10000:
    print(f"Turbulent flow: {Re_D:.2e} > 10000")
else:
    print(f"Transitional flow: 2300 < {Re_D:.2e} < 10000")


# Nusselt number and convection coefficient
# Ts < Tm (surface is colder than air), so n = 0.3 in the Dittus-Boelter equation
nu_avg_D  = nu_avg_turb(Re_D, Pr, 0.3)
h_avg = nu_avg_D * k / D    #W/m^2·K

print(f"Average Nusselt number: {nu_avg_D:.4f}")
print(f"Average convection coefficient: {h_avg:.4f} W/m^2·K")

# Solve for required tube length L
L = L_required(Ts, T_mi, T_mo, h_avg, P, m_dot, cp)
print(f"\nRequired tube length: L = {L:.2f} m")

# Verify L/D (length-to-diameter ratio) condition for the Dittus-Boelter equation
l_to_d = L / D
print(f"L/D = {l_to_d:.1f} {'(>= 10, Dittus-Boelter equation valid)' if l_to_d >= 10 else '(< 10, Dittus-Boelter equation may not be valid)'}")

# Time to melt ice
# Each tube extracts heat at rate q = m_dot * cp * (T_mi - T_mo)
q_per_tube = m_dot * cp * (T_mi - T_mo) #W, heat extracted per tube
q_total    = N_tubes * q_per_tube       #W, total heat extraction rate

print(f"\nHeat extracted per tube: {q_per_tube:.2f} W")
print(f"Total heat extraction rate: {q_total:.2f} W")

# Calculating mass, energy, and time of ice to melt
V_ice   = ice_frac * V_tank     #m^3, initial volume of ice
m_ice   = rho_ice * V_ice       #kg, mass of ice
E_melt  = m_ice * h_sf          #J, energy needed to melt all ice
t_melt  = E_melt / q_total      #s, time to melt ice
t_days  = t_melt / (3600 * 24)  #days, time to melt ice

print(f"\nVolume of ice: {V_ice:.1f} m^3")
print(f"Mass of ice: {m_ice:.1f} kg")
print(f"Energy to melt ice: {E_melt:.4e} J")
print(f"Time to melt ice: {t_melt:.1f} sec = {t_days:.1f} days")

# SIGNIFICANCE
# For air flowing through copper pipes submerged in an ice bath, the required length of the pipes is not that long (1.56 m) to cool the air by 10˚C, showing that turbulent internal convection is quite effective for cooling applications at relatively low mass flow rates.
