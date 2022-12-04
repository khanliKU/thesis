function [T1, m1] = tank_discharge_T_io(T0,m0,Tmd,md,gam)
%
% Calculate temperature after discharge
%
%   tank_discharge_T_io(T0,m0,md,gam)
%
%   T0: initial temperature in K
%   m0: initial gas mass in kg
%   Tmd: temperature of io gas in K
%   md: gas io in kg
%   gam: specific heat ratio Cp/Cv
    m1 = m0 + sum(md);
    T1 = (gam * md * Tmd' + m0 * T0)/m1;
end