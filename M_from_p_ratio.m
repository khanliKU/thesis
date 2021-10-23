function M = M_from_p_ratio(p_i,p_o,gam)
%
% Calculate flow velocity in Mach number
%
%   M_from_p_ratio(p_i,p_o,gam)
%
%   p_i: input pressure in Pa
%   p_o: output pressure in Pa
%   gam: specific heat ratio Cp/Cv
    M = sqrt(((p_i/p_o)^((gam-1)/gam)-1)*2/(gam-1));
end