function P = tank_discharge_P(P0,T0,m0,md,gam)
    T1 = tank_discharge_T(T0,m0,md,gam);
    P = P0 * T1 * (m0 - md) / (m0 * T0);
end