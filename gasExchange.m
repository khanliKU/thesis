function [rP1,rT1,rM1,rP2,rT2,rM2,mdot] = gasExchange(P1,T1,M1,P2,T2,M2,gam,A,R,dt,mdot_old)

    if P1 > P2
        mdot = checkImag(mdot_orifice(P1,P2,gam,A,T1,R));

        delta_m = 0.5*(mdot+mdot_old) * dt;
        
        if M1 == inf
            rT1 = T1;
            rM1 = M1;
        else
            [rT1, rM1] = tank_discharge_T_io(T1,M1,T1,-delta_m,gam);
        end
        if M2 == inf
            rT2 = T2;
            rM2 = M2;
        else
            [rT2, rM2] = tank_discharge_T_io(T2,M2,T1, delta_m,gam);
        end

        rP1 = P1 * (rT1/T1)^(gam/(gam-1));
        rP2 = P2 * (rT2/T2)^(gam/(gam-1));
    elseif P2 > P1
        mdot = -checkImag(mdot_orifice(P2,P1,gam,A,T2,R));

        delta_m = 0.5*(mdot+mdot_old) * dt;

        if M1 == inf
            rT1 = T1;
            rM1 = M1;
        else
            [rT1, rM1] = tank_discharge_T_io(T1,M1,T2,-delta_m,gam);
        end
        if M2 == inf
            rT2 = T2;
            rM2 = M2;
        else
            [rT2, rM2] = tank_discharge_T_io(T2,M2,T2, delta_m,gam);
        end

        rP1 = P1 * (rT1/T1)^(gam/(gam-1));
        rP2 = P2 * (rT2/T2)^(gam/(gam-1));
    else
        mdot = 0;
        rP1 = P1;
        rT1 = T1;
        rM1 = M1;
        rP2 = P2;
        rT2 = T2;
        rM2 = M2;
    end
    
end