function [Fmag,N,wire_len,wire_R,sol_V,sol_P,L,MR_total] = valve_magnetic_force(valve,gap,current)
    air_perm = 1.25663753*1e-6;
    %Magnetic Top
    magnetic_top_boss_center_cylnder_r = 0.5 * valve.magnetic_top_boss_dout;
    shell_center_cylinder_r = 0.25 * (valve.shell_dout + valve.shell_din);
    MR_volume_magnetic_top_hor_l = (shell_center_cylinder_r - magnetic_top_boss_center_cylnder_r) * 1e-3;
    MR_volume_magnetic_top_hor_A =...
        pi * (magnetic_top_boss_center_cylnder_r + shell_center_cylinder_r) * ... mean cylnder diameter
        valve.magnetic_top_h * 1e-6;
    MR_magnetic_top_hor = MR_volume_magnetic_top_hor_l / ( valve.perm * MR_volume_magnetic_top_hor_A );
    MR_volume_magnetic_top_ver_l_out = valve.magnetic_top_h * 1e-3;
    MR_volume_magnetic_top_ver_A_out = 0.25 * pi * (valve.magnetic_top_dout^2 - valve.shell_din^2) * 1e-6;
    MR_magnetic_top_ver_out = MR_volume_magnetic_top_ver_l_out / ( valve.perm * MR_volume_magnetic_top_ver_A_out ); % mm / (H/m mm^2) -> (m mm / H mm^2) -> 1e3 / H
    MR_volume_magnetic_top_ver_l_in = (MR_volume_magnetic_top_ver_l_out + valve.magnetic_top_boss_h) * 1e-3;
    MR_volume_magnetic_top_ver_A_in = 0.25 * pi * valve.magnetic_top_boss_dout^2 * 1e-6;
    MR_magnetic_top_ver_in = MR_volume_magnetic_top_ver_l_in / ( valve.perm * MR_volume_magnetic_top_ver_A_in ) ;% mm / (H/m mm^2) -> (m mm / H mm^2) -> 1e3 / H
    MR_magnetic_top = (MR_magnetic_top_hor + MR_magnetic_top_ver_out + MR_magnetic_top_ver_in);
    %Magnetic Bottom
    MR_surface_magnetic_bottom_hor_l = 0.5 * (valve.magnetic_bottom_dout - valve.magnetic_bottom_din) * 1e-3;
    MR_surface_magnetic_bottom_hor_A = pi * 0.5 * (valve.magnetic_bottom_dout + valve.magnetic_bottom_din) * valve.magnetic_bottom_h * 1e-6;
    MR_magnetic_bottom = MR_surface_magnetic_bottom_hor_l / (valve.perm * MR_surface_magnetic_bottom_hor_A);
    %Shell
    MR_surface_shell_hor_l = 0.5 * (valve.shell_dout - valve.shell_din) * 1e-3;
    MR_surface_shell_hor_A = pi * 0.5 * (valve.shell_dout + valve.shell_din) * valve.magnetic_bottom_h * 1e-6;
    MR_shell_hor = MR_surface_shell_hor_l / (valve.perm * MR_surface_shell_hor_A);
    MR_surface_shell_ver_l = (valve.shell_h - 0.5 * valve.magnetic_bottom_h) * 1e-3;
    MR_surface_shell_ver_A = pi * 0.25 * (valve.shell_dout^2 -valve.shell_din^2) * 1e-6;
    MR_shell_ver = MR_surface_shell_ver_l / (valve.perm * MR_surface_shell_ver_A);
    MR_shell = MR_shell_ver + MR_shell_hor;
    %Shell - Magnetic Bottom clearance
    MR_surface_shell_mag_bot_hor_l = 0.5 * valve.clearance * 1e-3;
    MR_surface_shell_mag_bot_hor_A = pi * 0.5 * (valve.shell_din + valve.magnetic_bottom_dout) * valve.magnetic_bottom_h * 1e-6;
    MR_shell_mag_bot_hor = MR_surface_shell_mag_bot_hor_l / (valve.perm * MR_surface_shell_mag_bot_hor_A);
    %Shell - Magnetic Top clearance
    %MR_surface_shell_mag_top_hor_l = 0.5 * valve.clearance * 1e-3
    %MR_surface_shell_mag_top_hor_A = pi * 0.5 * (valve.shell_din + valve.magnetic_top_dout)...
    %    * (valve.magnetic_top_h - valve.magnetic_top_hole_h) * 1e-6
    MR_shell_mag_top_hor = 0; %MR_surface_shell_mag_top_hor_l / (valve.perm * MR_surface_shell_mag_top_hor_A)
    %Gap
    MR_gap_hor_l = gap * 1e-3;
    MR_gap_hor_A = 0.25 * pi * valve.magnetic_core_d^2 * 1e-6;
    MR_gap = MR_gap_hor_l / (air_perm * MR_gap_hor_A);
    %Valve Spool
    MR_surface_valve_spool_hor_l = valve.magnetic_core_d * 1e-3;
    MR_surface_valve_spool_hor_A = pi * valve.magnetic_core_d * valve.magnetic_bottom_h * 1e-6;
    MR_valve_spool_hor = MR_surface_valve_spool_hor_l / (valve.perm * MR_surface_valve_spool_hor_A);
    MR_surface_valve_spool_ver_l = (0.5 * valve.shell_h - gap) * 1e-3;
    MR_surface_valve_spool_ver_A = pi * 0.25 * valve.magnetic_core_d^2 * 1e-6;
    MR_valve_spool_ver = MR_surface_valve_spool_ver_l / (valve.perm * MR_surface_valve_spool_ver_A);
    MR_valve_spool = MR_valve_spool_ver + MR_valve_spool_hor;
    %Valve Spool - Magnetic Bottom valve.clearance
    MR_surface_valve_spool_magnetic_bottom_hor_l = 0.5 * valve.clearance * 1e-3;
    MR_surface_valve_spool_magnetic_bottom_hor_A = pi * 0.5 * (valve.magnetic_core_d + valve.magnetic_bottom_din) * valve.magnetic_bottom_h * 1e-6;
    MR_valve_spool_magnetic_bottom_hor = MR_surface_valve_spool_magnetic_bottom_hor_l / (air_perm * MR_surface_valve_spool_magnetic_bottom_hor_A);
    %Total
    MR_total = MR_valve_spool + MR_gap + MR_shell + MR_magnetic_bottom + MR_magnetic_top +...
        MR_valve_spool_magnetic_bottom_hor + MR_shell_mag_bot_hor + MR_shell_mag_top_hor;
    R_per_km = 18.426905*valve.wire_area^-0.997135;
    wire_d = 2 * sqrt(valve.wire_area/pi);
    coil_cross_section_A = 0.5 * (valve.coil_dout - valve.coil_din) * valve.coil_h; % mm2;
    N = coil_cross_section_A / wire_d^2; % TODO remove safety factor
    wire_len = N * pi * (valve.coil_dout + valve.coil_din) * 1e-3;
    wire_R = wire_len * R_per_km * 1e-3;
    sol_V = wire_R * current;
    sol_P = wire_R * current^2;
    flux = N * current / MR_total;
    Fmag = 0.5 * flux^2 / (air_perm * MR_gap_hor_A);
    % L = N * flux / current;
    L = N^2 / MR_total; % This derivation lets me compute inductance in terms of air gap, not current.
    %figure
    %plot([MR_valve_spool,MR_gap,MR_shell,MR_magnetic_bottom,MR_magnetic_top,...
    %    MR_valve_spool_magnetic_bottom_hor,MR_shell_mag_bot_hor,MR_shell_mag_top_hor])
    % valve_magnetic_force(dynamic.valve,0.023,1)
end