function [d, F] = dFmag_dd4(l1,l2,l3,l4,lc,d1,d2,d3,d4)
    mu0 = 1.25663706212E-6;
    mu = 850 * mu0;
    la = l1/2 + l3/2 + l2;

    packingFactor = 0.85;
    A_wire = 0.22E-6;

    rTot = theoricTotalReluctance(l1,l2,l3,lc,d1,d2,d4)

    S2 = pi * d4 /4 ;
    N = floor(l4 * (d2 - d3) * packingFactor / A_wire)
    l_wire = N .* (d2 + d3) .* pi
    R_wire = l_wire * 39E-3
    I = 2;

    flux = N * I /rTot;

    d = 1/(2*mu0)...
        *(2 * (...
        -(l4 * packingFactor / A_wire) *...
        rTot - N .* I .*...
        (-2 *...
        (mu0 * (la-lc) + mu * lc) ./ (mu * mu0 * 0.25 * pi * d4.^3)))...
        /(rTot.^2) *...
        flux * 1./S2 - 0.5 * flux^2 * pi * d4 * 1/S2.^2);

    F = 0.5 * flux.^2 ./ (mu0 * S2);
end