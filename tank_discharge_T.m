function T1 = tank_discharge_T(T0,m0,md,gam)
%
% Calculate temperature after discharge
%
%   tank_discharge_T(T0,m0,md,gam)
%
%   T0: initial temperature in K
%   m0: initial gas mass in kg
%   md: gas release in kg
%   gam: specific heat ratio Cp/Cv
    T1 = (m0 * T0 - 0.5 * md * gam * T0)/...
        (m0 - md + 0.5 * md * gam);
end