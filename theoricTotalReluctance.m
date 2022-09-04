function R = theoricTotalReluctance(l1,l2,l3,lc,d1,d2,d4)
    mu0 = 1.25663706212E-6;
    mu = 850 * mu0;
    la = l1/2 + l3/2 + l2;
    lb = 0.5 * (d1 + d2);

    S1 = 0.25*pi*(d1^2-d2^2);
    S2 = 0.25*pi*d4.^2;
    S3 = 0.25*pi*l1*(d1+d2);

    Rvert_in = (la-lc)./(mu*S2) + lc./(mu0*S2);
    Rvert_out = la/(mu*S1);
    Rhor_top = lb/(mu*S3);

    R = Rvert_out + Rvert_in + 2 * Rhor_top;
end