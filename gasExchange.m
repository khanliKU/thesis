function [rP1,rT1,rM1,rP2,rT2,rM2,mdot] = gasExchange(P1,T1,M1,vol1,P2,T2,M2,vol2,gam,A,R,Cv,dt,mdot_old)

    mdot = 0;
    vel = 0;
    rP1 = P1;
    rT1 = T1;
    rM1 = M1;
    rP2 = P2;
    rT2 = T2;
    rM2 = M2;

    if A == 0
        return
    end

    if P1 > P2
        rho = M1/vol1;
        if isnan(rho)
            rho = P1 / (T1 * R);
        end
        vel = M_orifice(P1,P2,gam) * speed_of_sound(gam,P1,rho);
        mdot = checkImag(mdot_orifice(P1,P2,gam,A,T1,R));
    elseif P2 > P1
        rho = M2/vol2;
        if isnan(rho)
            rho = P2 / (T2 * R);
        end
        vel = M_orifice(P2,P1,gam) * speed_of_sound(gam,P2,rho);
        mdot = -checkImag(mdot_orifice(P2,P1,gam,A,T2,R));
    end

    mdot_ave = mdot;

    if mdot_ave > 0
        % Flow from Vol 1 to Vol 2
        delta_m = mdot_ave * dt;
        
        if M1 == inf
            rP1 = P1;
            rT1 = T1;
            rM1 = M1;
        else
            [rP1, rT1, rM1] = tank_discharge_io_pv(T1,M1,P1,vel*A,-delta_m,R,vol1,Cv);
        end
        if M2 == inf
            rP2 = P2;
            rT2 = T2;
            rM2 = M2;
        else
            [rP2, rT2, rM2] = tank_discharge_io_pv(T2,M2,P1,vel*A, delta_m,R,vol2,Cv);
        end

    elseif mdot_ave < 0
        % Flow from Vol 2 to Vol 1
        delta_m = mdot_ave * dt;
        
        if M1 == inf
            rP1 = P1;
            rT1 = T1;
            rM1 = M1;
        else
            [rP1, rT1, rM1] = tank_discharge_io_pv(T1,M1,P2,vel*A,-delta_m,R,vol1,Cv);
        end
        if M2 == inf
            rP2 = P2;
            rT2 = T2;
            rM2 = M2;
        else
            [rP2, rT2, rM2] = tank_discharge_io_pv(T2,M2,P2,vel*A, delta_m,R,vol2,Cv);
        end
    end
    
    checkImag(mdot);
    checkImag(vel);
    checkImag(rP1);
    checkImag(rT1);
    checkImag(rM1);
    checkImag(rP2);
    checkImag(rT2);
    checkImag(rM2);
    if rP1 < 0 || rP2 < 0
        asda =12312;
    end
end