function [P, T, m1] = tank_discharge_io(P0,T0,m0,Tother,Pother,orifices,gam,R,dt)
%
% Calculate pressure after discharge
%
%   tank_discharge_P_io(P0,T0,m0,md,gam)
%
%   P0: initial pressure in Pa
%   T0: initial temperature in K
%   m0: initial gas mass in kg
%   Tmd: temperature of io gas in K
%   md: gas io in kg
%   gam: specific heat ratio Cp/Cv
    md = mdot_orifice(Pother,P0,gam,orifices,Tother,R)*dt; % kg/s;
    [T, m1] = tank_discharge_T_io(T0,m0,Tother,md,gam);
    P = P0 * T * m1 / (m0 * T0);
end