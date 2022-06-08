function [F, P, I, N, ra, ro, wLen, Wf, I_safe] = coil_model(arm_perm, r0, s, l, t, V, x, Aw, Rwm)
% Inputs
% arm_perm: relative permeability of armature
% c0: coil inner diameter
% s: coil turn stack height
% l: coil length
% t: winding type l or g
% V: voltage applied to coil
% x: armature displacement
% Ouputs
% F: Force in N
% P: Power in W
% I: Current in Amps
% N: # of turns
% ra: mean radius
% ro: outer radius
% wLen: wire length
% Wf: winding factor
% I_safe: Current in Amps
    rec_max_A = 3.5; % A/mm2
    [Wf, ra, ro, wLen, N] = winding_from_r0_s_l_Aw_t(r0,s,l,Aw,t);
    A = pi*ra^2*1e-6;
    L = coil_inductance(arm_perm,N,A,l,x);
    R = wLen * Rwm;
    I = V/R;
    I_safe = Aw * rec_max_A;
    F = coil_force(V,R,L,log(arm_perm),l);
    P = V*I;
end