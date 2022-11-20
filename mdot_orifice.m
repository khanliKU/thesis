function mdot = mdot_orifice(p_i,p_o,gam,A,Tt,R)
%
% Calculate mass flow rate through orifice
%
%   mdot_orifice(p_i,p_o,gam,A,Tt,R)
%
%   p_i: input pressure in Bar
%   p_o: output pressure in Bar
%   gam: specific heat ratio Cp/Cv
%   A: orifice area in m2
%   Tt: initial temperature in K
%   R: individual gas constant J/KgK
    M = M_orifice(p_i,p_o,gam);
    mdot = A * B_to_Pa(p_i) / sqrt(Tt) * sqrt(gam/R) * M *...
        (1 + (gam-1)/2 * M^2)^...
        -((gam+1)/(2*(gam-1)));
end