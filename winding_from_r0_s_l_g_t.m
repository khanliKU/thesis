function [Wf, ra, ro, wLen, N] = winding_from_r0_s_l_Aw_t(r0,s,l,Aw,t)
% Inputs
% s: stack count
% r0: innter r
% l: coil length
% Aw: wire cross section area mm^2
% t: type
% Ouputs
% Wf: winding factor
% ra: mean radius
% ro: outer radius
% wLen: wire length
% N: # of turns
    dw = sqrt(4*Aw/pi);
    if s > 1
        if t == 'l' % lattice
            ro = r0 + s * dw * 1/sqrt(3);
        elseif t == 'g' % grid
            ro = r0 + s * dw;
        end
        ra = 0.5 * (r0 + ro);
        N = s * l / dw ;
        wLen = 2 * pi * ra * l / dw;
        Wf = (r0/ra)^2;
    end
end