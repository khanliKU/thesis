%% Inline Double Plunger 00

clearvars
%% Outlet

goal_Cv = 6
goal_Cd = 0.61
goal_B = 0
goal_r = 0.00465*sqrt(goal_Cv/goal_Cd) %m
goal_A = pi * goal_r^2 %m2
Cd_cylinder = 0.5
%% Main

main_orifice_r = 2e-3 %m
main_tower_t = 2e-3 %m
main_tower_r = main_orifice_r + main_tower_t %m
main_tower_A = pi*main_tower_r^2 %m2
main_tower_h = 30e-3 %m
%% Outlet

outlet_orifice_A = goal_A*sqrt(goal_Cd/Cd_cylinder) + main_tower_A %m2
outlet_orifice_r = sqrt(outlet_orifice_A/pi) %m
outlet_orifice_h = 10e-3 %m
outlet_entrance_A = outlet_orifice_A * 4 %m2
outlet_entrance_r = sqrt(outlet_entrance_A/pi) %m
outlet_entrance_h = 10e-3 %m
outlet_boss_t = 5e-3 %m
outlet_boss_r = outlet_entrance_r + outlet_boss_t %m
outlet_boss_A = pi*outlet_boss_r^2 %m2
outlet_boss_h = outlet_orifice_h + outlet_entrance_h %m
outlet_t = 10e-3 %m
%% Main

main_A = outlet_boss_A + outlet_entrance_A %m2
main_r = sqrt(main_A/pi) %m
main_high_pressure_A = pi * (main_r^2 - outlet_boss_r^2) %m2
main_down_side_A = outlet_entrance_A %m2
main_low_pressure_A = (main_high_pressure_A + main_down_side_A)*1.00 %m2
main_boss_out_r = sqrt(main_low_pressure_A/pi) %m
main_boss_t = 3e-3 %m
main_boss_inner_r = main_boss_out_r - main_boss_t %m
% Main Force Balance

%{
HP = 1.5 %MPa
LP = 1.5 %MPA
% Up direction is positive
F_hp = HP * main_high_pressure_A * 1e3 %kN
F_lp = LP * main_low_pressure_A * 1e3 %kN
F_ds = LP * main_down_side_A * 1e3 %kN
F = F_hp - F_lp + F_ds %kN
% Pilot Open

step = 0.1
[UP,DP] = meshgrid(0.1:step:40,0:1:1);
figure(1)
F = (UP * main_high_pressure_A + UP .* DP * (main_down_side_A - main_low_pressure_A)) * 1e3; %kN
min(min(F))
contour(UP,DP,F,'ShowText','on')
xlabel('UP'),ylabel('DP')
grid on, grid minor
title('Pilot open')

figure(2)
mesh(UP,DP,F)
xlabel('UP [MPa]'),ylabel('DP %'),zlabel('F [kN]')
title('Pilot open')
% Pilot Closed

figure(3)
F = (UP * (main_high_pressure_A - main_low_pressure_A) + UP .* DP * main_down_side_A) * 1e3; %kN
min(min(F))
contour(UP,DP,F,'ShowText','on')
xlabel('UP'),ylabel('DP')
grid on, grid minor
title('Pilot closed')

figure(4)
mesh(UP,DP,F)
xlabel('UP [MPa]'),ylabel('DP %'),zlabel('F [kN]')
title('Pilot closed')
%}
%% Main

main_t = 10e-3 %m
main_movement_range = outlet_entrance_A / (2 * pi * outlet_entrance_r) %m
%% Retainer

retainer_main_groove_depth = main_t + main_movement_range + 2e-3 %m
retainer_main_groove_r = main_r %m
retainer_t = 4e-3 %m
retainer_r = retainer_main_groove_r + retainer_t %m
retainer_main_boss_groove_r = main_boss_out_r %m
%% Pilot

pilot_seal_r = main_boss_inner_r * 0.5 %m
pilot_seal_t = 4e-3 %m
pilot_stem_r = pilot_seal_r * 0.25 %m
pilot_steam_A = pilot_stem_r^2 * pi %m2
pilot_seal_A = pilot_seal_r^2 * pi %m2
pilot_movement_range = 20e-3 %m
%% Main

main_boss_h = pilot_movement_range * 1.05 %m
%% Retainer

retainer_main_boss_groove_depth = main_movement_range + 2e-3 %m
retainer_total_h = retainer_main_boss_groove_depth + retainer_main_groove_depth + retainer_t %m
%% Inlet

inlet_retainer_groove_r = sqrt(outlet_entrance_r^2 + retainer_r^2) %m
inlet_t = 8e-3 %m
inlet_r = inlet_retainer_groove_r + inlet_t %m
inlet_outlet_distance = outlet_entrance_A / (2 * pi * retainer_main_boss_groove_r) %m
inlet_h = 2 * inlet_outlet_distance + retainer_total_h + inlet_t %m