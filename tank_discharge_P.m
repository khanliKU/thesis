function [P, T] = tank_discharge_P(P0,T0,m0,md,gam)
%
% Calculate pressure after discharge
%
%   tank_discharge_P(P0,T0,m0,md,gam)
%
%   P0: initial pressure in Pa
%   T0: initial temperature in K
%   m0: initial gas mass in kg
%   md: gas release in kg
%   gam: specific heat ratio Cp/Cv
    T = tank_discharge_T(T0,m0,md,gam);
    P = P0 * T * (m0 - md) / (m0 * T0);
end