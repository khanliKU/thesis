function F = coil_force(V,R,L,a,l)
    F = -0.5*(V/R)^2*a/l*L;
end