function M = M_orifice(p_i,p_o,gam)
    M = M_from_p_ratio(p_i,p_o,gam);
    if M > 1
        M = 1;
    end
end