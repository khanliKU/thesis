function mdot = mdot_orifice(p_i,p_o,gam,A,T_i,T_o,R)
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
    for i = 1:length(p_i)
        if p_i(i) > p_o
            p_h(i) = p_i(i);
            p_l(i) = p_o;
            Tt(i) = T_i(i);
        else
            p_h(i) = p_o;
            p_l(i) = p_i(i);
            Tt(i) = T_o;
        end
    end
    M = M_orifice(p_h,p_l,gam);
    mdot = A .* B_to_Pa(p_i) ./ sqrt(Tt) .* sqrt(gam/R) .* M .*...
        (1 + (gam-1)/2 * M.^2).^...
        -((gam+1)/(2*(gam-1)));
    for i = 1:length(p_i)
        if p_i(i) < p_o
            mdot(i) = - mdot(i);
        elseif p_i(i) == p_o
            mdot(i) = 0;
        end
    end
end