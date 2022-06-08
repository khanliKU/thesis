function L = coil_inductance(arm_perm,N,A,l,x)
    air_perm = 4*pi*1e-7; % H/m
    L0 = arm_perm * air_perm * N^2 * A / (l * 1e-3);
    L = L0*exp(-log(arm_perm)*x/l);
end