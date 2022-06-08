import numpy as np
#
# Calculate pressure after discharge
#
#   tank_discharge_P(P0,T0,m0,md,gam)
#
#   P0: initial pressure in Pa
#   T0: initial temperature in K
#   m0: initial gas mass in kg
#   md: gas release in kg
#   gam: specific heat ratio Cp/Cv
def tank_discharge_P(P0,T0,m0,md,gam):
    T1 = tank_discharge_T(T0,m0,md,gam)
    P = P0 * T1 * (m0 - md) / (m0 * T0)
    return P, T1

#
# Calculate temperature after discharge
#
#   tank_discharge_T(T0,m0,md,gam)
#
#   T0: initial temperature in K
#   m0: initial gas mass in kg
#   md: gas release in kg
#   gam: specific heat ratio Cp/Cv
def tank_discharge_T(T0,m0,md,gam):
    T1 = (m0 * T0 - 0.5 * md * gam * T0)\
        (m0 - md + 0.5 * md * gam)
    return T1

#
# Calculate mass flow rate through orifice
#
#   mdot_orifice(p_i,p_o,gam,A,Tt,R)
#
#   p_i: input pressure in Pa
#   p_o: output pressure in Pa
#   gam: specific heat ratio Cp/Cv
#   A: orifice area in m2
#   Tt: initial temperature in K
#   R: individual gas constant J/KgK
def mdot_orifice(p_i,p_o,gam,A,Tt,R):
    M = M_orifice(p_i,p_o,gam)
    mdot = A * B_to_Pa(p_i) / np.sqrt(Tt) * np.sqrt(gam/R) * M *\
        (1 + (gam+1)/2 * M**2)**\
        -((gam+1)/(2*(gam-1)))
    return mdot

#
# Calculate flow velocity in Mach number
#
#   M_from_p_ratio(p_i,p_o,gam)
#
#   p_i: input pressure in Pa
#   p_o: output pressure in Pa
#   gam: specific heat ratio Cp/Cv
def M_from_p_ratio(p_i,p_o,gam):
    M = np.sqrt(((p_i/p_o)**((gam-1)/gam)-1)*2/(gam-1))
    return M

#
# Calculate flow velocity in Mach number through orifice
#
#   M_orifice(p_i,p_o,gam)
#
#   p_i: input pressure in Pa
#   p_o: output pressure in Pa
#   gam: specific heat ratio Cp/Cv
def M_orifice(p_i,p_o,gam):
    M = M_from_p_ratio(p_i,p_o,gam)
    if M > 1:
        M = 1
    return M

#
# Convert Bar to Pascal
#
#   B_to_Pa(B)
#
#   B: pressure in bar
def B_to_Pa(B):
    Pa = 101300 * B
    return Pa

#
# Calculate mass from ideal gas equation
#
#   ideal_gas_mass(P,V,T,r)
#
#   P: pressure in Pa
#   V: volume in m3
#   T: temperature in K
#   r: individual gas constant J/KgK
def ideal_gas_mass(P,V,T,r):
    m = P * V / (r * T)
    return m

class Gas:
    def __init__(self, gamma, R):
        self.gamma = gamma
        self.R = R

class Volume:
    def __init__(self, gas, volume, pressure, temperature, mass):
        self.mass = mass
        self.volume = volume
        self.pressure = pressure
        self.temperature = temperature
        self.mass = mass
    
    #
    # Calculate state of the volume in next iteration
    #
    #   io(delta_m)
    #
    #   delta_m: change in mass in kg
    def io(self,delta_m):
        Pnext, Tnext = tank_discharge_P(self.pressure,self.temperature,\
                                        self.mass,delta_m,self.gas.gamma)
        return Volume(self.gas,self.volume,Pnext,Tnext,self.mass+delta_m)


































