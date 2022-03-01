function T1 = tank_discharge_T_io(T0,m0,Tmd,md,gam)
%
% Calculate temperature after discharge
%
%   tank_discharge_T(T0,m0,md,gam)
%
%   T0: initial temperature in K
%   m0: initial gas mass in kg
%   Tmd: temperature of io gas in K
%   md: gas io in kg
%   gam: specific heat ratio Cp/Cv
    T1 = (md * gam * Tmd + m0 * T0)/...
        (m0 + md);
end