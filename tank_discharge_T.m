function T1 = tank_discharge_T(T0,m0,md,gam)
    T1 = (m0 * T0 - 0.5 * md * gam * T0)/...
        (m0 - md + 0.5 * md * gam);
end