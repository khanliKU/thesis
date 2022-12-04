function F = force_micro_hs_v(rel_perm, gap, N, I)
%https://www.electrical4u.com/magnetic-reluctance/
filename = 'design_params_ref.xlsx';
[NUM,TXT,RAW] = xlsread(filename);
PI = pi;
for i=1:size(RAW,1)
    eval(RAW(i,1) + " = " + RAW(i,2))
end

%Sealin Surface
do = 0.5;
r = 0.25;
R = 0.5 * do + r;
fun = @(x) 2*pi*(sqrt(r^2 - x.^2) + R).*sqrt(r^2./(r^2 - x.^2))
S = integral(fun,0,r)
Ao = pi * do^2 / 4
Ap = pi * R^2
20.5*S
50*Ap

%Magnetic Reluctance
%rel_perm = 850
air_perm = 1.25663753*1e-6
perm = rel_perm * air_perm
%Magnetic Top
magnetic_top_boss_center_cylnder_r = 0.25 * (magnetic_top_boss_dout + magnetic_top_boss_din)
shell_center_cylinder_r = 0.25 * (shell_dout + shell_din)
MR_volume_magnetic_top_hor_l = (shell_center_cylinder_r - magnetic_top_boss_center_cylnder_r) * 1e-3
MR_volume_magnetic_top_hor_A =...
    (pi * (magnetic_top_boss_center_cylnder_r + shell_center_cylinder_r) * ... mean cylnder diameter
    (magnetic_top_h - magnetic_top_hole_h)) * 1e-6
MR_magnetic_top_hor = MR_volume_magnetic_top_hor_l / ( perm * MR_volume_magnetic_top_hor_A )
MR_volume_magnetic_top_ver_l_out = 0.5 * (magnetic_top_h - magnetic_top_hole_h) * 1e-3
MR_volume_magnetic_top_ver_A_out = 0.25 * pi * (magnetic_top_dout^2 - shell_din^2) * 1e-6
MR_magnetic_top_ver_out = MR_volume_magnetic_top_ver_l_out / ( perm * MR_volume_magnetic_top_ver_A_out ) % mm / (H/m mm^2) -> (m mm / H mm^2) -> 1e3 / H
MR_volume_magnetic_top_ver_l_in = (MR_volume_magnetic_top_ver_l_out + magnetic_top_boss_h) * 1e-3
MR_volume_magnetic_top_ver_A_in = 0.25 * pi * (magnetic_top_boss_dout^2 - magnetic_top_boss_din^2) * 1e-6
MR_magnetic_top_ver_in = MR_volume_magnetic_top_ver_l_in / ( perm * MR_volume_magnetic_top_ver_A_in ) % mm / (H/m mm^2) -> (m mm / H mm^2) -> 1e3 / H
MR_magnetic_top = (MR_magnetic_top_hor + MR_magnetic_top_ver_out + MR_magnetic_top_ver_in)
mmf = (200*0.75)*MR_magnetic_top
%Magnetic Bottom
MR_surface_magnetic_bottom_hor_l = 0.5 * (magnetic_bottom_dout - magnetic_bottom_din) * 1e-3
MR_surface_magnetic_bottom_hor_A = pi * 0.5 * (magnetic_bottom_dout + magnetic_bottom_din) * magnetic_bottom_h * 1e-6
MR_magnetic_bottom = MR_surface_magnetic_bottom_hor_l / (perm * MR_surface_magnetic_bottom_hor_A)
%Shell
MR_surface_shell_hor_l = 0.5 * (shell_dout - shell_din) * 1e-3
MR_surface_shell_hor_A = pi * 0.5 * (shell_dout + shell_din) * magnetic_bottom_h * 1e-6
MR_shell_hor = MR_surface_shell_hor_l / (perm * MR_surface_shell_hor_A)
MR_surface_shell_ver_l = (shell_h - 0.5 * magnetic_bottom_h) * 1e-3
MR_surface_shell_ver_A = pi * 0.25 * (shell_dout^2 -shell_din^2) * 1e-6
MR_shell_ver = MR_surface_shell_ver_l / (perm * MR_surface_shell_ver_A)
MR_shell = MR_shell_ver + MR_shell_hor
%Gap
%gap = 0.2
MR_gap_hor_l = gap * 1e-3
MR_gap_hor_A = 0.25 * pi * valve_spool_d^2 * 1e-6
MR_gap = MR_gap_hor_l / (air_perm * MR_gap_hor_A)
%Valve Spool
MR_surface_valve_spool_hor_l = 0.5 * (valve_spool_d - valve_spool_hole_d) * 1e-3
MR_surface_valve_spool_hor_A = pi * 0.5 * (valve_spool_d + valve_spool_hole_d) * magnetic_bottom_h * 1e-6
MR_valve_spool_hor = MR_surface_valve_spool_hor_l / (perm * MR_surface_valve_spool_hor_A)
MR_surface_valve_spool_ver_l = (0.5 * shell_h - gap) * 1e-3
MR_surface_valve_spool_ver_A = pi * 0.25 * (valve_spool_d^2 - valve_spool_hole_d^2) * 1e-6
MR_valve_spool_ver = MR_surface_valve_spool_ver_l / (perm * MR_surface_valve_spool_ver_A)
MR_valve_spool = MR_valve_spool_ver + MR_valve_spool_hor
%Total
MR_total = MR_valve_spool + MR_gap + MR_shell + MR_magnetic_bottom + MR_magnetic_top
NI = N * I % N: # of turns, I: current per turn
flux = NI / MR_total
F = flux^2 / (air_perm * MR_gap_hor_A)
end
