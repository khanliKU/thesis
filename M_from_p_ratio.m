function M = M_from_p_ratio(p_i,p_o,gam)
    M = sqrt(((p_i/p_o)^((gam-1)/gam)-1)*2/(gam-1));
end