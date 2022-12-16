function [P1, T1, m1] = tank_discharge_io_pv(T0,m0,P_md,v_md,md,R,vol,Cv)
%
% Calculate temperature after discharge
%
%   tank_discharge_T_io(T0,m0,md,gam)
%
%   T0: initial temperature in K
%   m0: initial gas mass in kg
%   P_md: Pressure of io gas in Pa
%   v_md: specific volume of io gas in kg/m3
%   md: gas io in kg
%   R: Gas constant in kJ/kg.K
%   vol: volume of tank
%   Cv: specific heat constant volume in kJ/kg.K
    m1 = m0 + sum(md);
    if m1 < 0
        a = 13;
    end
    h = P_md .* v_md .* md;
    u0 = m0 * Cv * T0;
    u1 = u0 + sum(h);
    T1 = u1 / (m1 * Cv);
    P1 = R * T1 / (vol/m1);
end