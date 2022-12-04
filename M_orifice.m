function M = M_orifice(p_i,p_o,gam)
%
% Calculate flow velocity in Mach number through orifice
%
%   M_orifice(p_i,p_o,gam)
%
%   p_i: input pressure in Pa
%   p_o: output pressure in Pa
%   gam: specific heat ratio Cp/Cv
    M = M_from_p_ratio(p_i,p_o,gam);
    if M > 1
        M = 1;
    end
end