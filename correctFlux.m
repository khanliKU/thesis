function r = correctFlux(magFieldStr_in_Am,bh430F)
    magFieldStr = magFieldStr_in_Am * (4*pi) / 1e3;
    if magFieldStr > bh430F.H_Oe(end)
        r_in_gauss = (bh430F.B_Gauss(end)-bh430F.B_Gauss(end-1))*(magFieldStr-bh430F.H_Oe(end))+bh430F.B_Gauss(end);
    else
        r_in_gauss = interp1(bh430F.H_Oe,bh430F.B_Gauss,magFieldStr);
    end
    r = r_in_gauss / 1e4;
end