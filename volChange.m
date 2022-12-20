function [rT, rP] = volChange(V1,V2,T,P,gam)
    if V1 == V2
        rT = T;
        rP = P;
    else
        rT = checkImag(T * (V1/V2)^(gam-1));
        rP = checkImag(P * (V1/V2)^gam);
    end
end