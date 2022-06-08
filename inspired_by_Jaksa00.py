import numpy as np
import matplotlib.pyplot as plt
## Inline Double Plunger 00

## Outlet

goal_Cv = 1
goal_Cd = 0.61
goal_B = 0
goal_r = 0.00465*np.sqrt(goal_Cv/goal_Cd) #m
goal_A = np.pi * goal_r**2 #m2
Cd_cylinder = 0.5
## Main

peek_orifice_r = 1e-3 #m
peek_tower_t = 3e-3 #m
peek_tower_r = peek_orifice_r + peek_tower_t #m
peek_tower_A = np.pi*peek_tower_r**2 #m2
peek_tower_h = 10e-3 #m
## Outlet

outlet_orifice_A = goal_A*np.sqrt(goal_Cd/Cd_cylinder) + peek_tower_A #m2
outlet_orifice_r = np.sqrt(outlet_orifice_A/np.pi) #m
outlet_orifice_h = 10e-3 #m
outlet_entrance_A = outlet_orifice_A * 4 #m2
outlet_entrance_r = np.sqrt(outlet_entrance_A/np.pi) #m
outlet_entrance_h = 10e-3 #m
outlet_boss_t = 5e-3 #m
outlet_boss_r = outlet_entrance_r + outlet_boss_t #m
outlet_boss_A = np.pi*outlet_boss_r**2 #m2
outlet_boss_h = outlet_orifice_h + outlet_entrance_h #m
outlet_t = 10e-3 #m
## Main

main_fully_open_min_h = goal_A*np.sqrt(goal_Cd/Cd_cylinder) / (2*outlet_entrance_r)
main_working_volume = main_fully_open_min_h * outlet_entrance_A
main_A = outlet_boss_A + outlet_entrance_A #m2
main_r = np.sqrt(main_A/np.pi) #m
main_high_pressure_A = np.pi * (main_r**2 - outlet_boss_r**2) #m2
main_down_side_A = outlet_entrance_A #m2
main_low_pressure_A = (main_high_pressure_A + main_down_side_A)*1.00 #m2
main_boss_out_r = np.sqrt(main_low_pressure_A/np.pi) #m
main_boss_t = 3e-3 #m
main_boss_inner_r = main_boss_out_r - main_boss_t #m
# Main Force Balance

outlet_A_arr = []
main_fully_open_min_h_arr = []
main_working_volume_arr = []
for a in np.arange(outlet_orifice_A,outlet_entrance_A, (outlet_entrance_A - outlet_orifice_A ) / 2 ):
    outlet_A_arr.append(a)
    main_fully_open_min_h_arr.append(goal_A*np.sqrt(goal_Cd/Cd_cylinder) / (2*np.sqrt(a)/np.pi))
    main_working_volume_arr.append(main_fully_open_min_h_arr[-1] * a)
    
plt.figure(0)
plt.plot(outlet_A_arr,main_working_volume_arr)
plt.figure(1)
plt.plot(outlet_A_arr,main_fully_open_min_h_arr)
#{
HP = 1.5 #MPa
LP = 1.5 #MPA
# Up direction is positive
F_hp = HP * main_high_pressure_A * 1e3 #kN
F_lp = LP * main_low_pressure_A * 1e3 #kN
F_ds = LP * main_down_side_A * 1e3 #kN
F = F_hp - F_lp + F_ds #kN
# pilot Open

# step = 0.1
# #[UP,DP] = meshgrid(0.1:step:40,0:1:1);
# figure(1)
# F = (UP * main_high_pressure_A + UP .* DP * (main_down_side_A - main_low_pressure_A)) * 1e3; #kN
# min(min(F))
# contour(UP,DP,F,'ShowText','on')
# xlabel('UP'),ylabel('DP')
# grid on, grid minor
# title('pilot open')

# figure(2)
# mesh(UP,DP,F)
# xlabel('UP [MPa]'),ylabel('DP #'),zlabel('F [kN]')
# title('pilot open')
# # pilot Closed

# figure(3)
# F = (UP * (main_high_pressure_A - main_low_pressure_A) + UP .* DP * main_down_side_A) * 1e3; #kN
# min(min(F))
# contour(UP,DP,F,'ShowText','on')
# xlabel('UP'),ylabel('DP')
# grid on, grid minor
# title('pilot closed')

# figure(4)
# mesh(UP,DP,F)
# xlabel('UP [MPa]'),ylabel('DP #'),zlabel('F [kN]')
# title('pilot closed')
# #}
## Main

main_t = 10e-3 #m
main_movement_range = outlet_entrance_A / (2 * np.pi * outlet_entrance_r) #m
## Retainer

retainer_main_groove_depth = main_t + main_movement_range + 2e-3 #m
retainer_main_groove_r = main_r #m
retainer_t = 4e-3 #m
retainer_r = retainer_main_groove_r + retainer_t #m
retainer_main_boss_groove_r = main_boss_out_r #m
## pilot

pilot_seal_r = main_boss_inner_r * 0.5 #m
pilot_seal_t = 4e-3 #m
pilot_stem_r = pilot_seal_r * 0.25 #m
pilot_stem_A = pilot_stem_r**2 * np.pi #m2
pilot_seal_A = pilot_seal_r**2 * np.pi #m2
pilot_movement_range = 20e-3 #m
## Main

main_boss_h = pilot_movement_range * 1.05 #m
## Retainer

retainer_main_boss_groove_depth = main_movement_range + 2e-3 #m
retainer_total_h = retainer_main_boss_groove_depth + retainer_main_groove_depth + retainer_t #m
## Inlet

inlet_retainer_groove_r = np.sqrt(outlet_entrance_r**2 + retainer_r**2) #m
inlet_t = 8e-3 #m
inlet_r = inlet_retainer_groove_r + inlet_t #m
inlet_outlet_distance = outlet_entrance_A / (2 * np.pi * retainer_main_boss_groove_r) #m
inlet_h = 2 * inlet_outlet_distance + retainer_total_h + inlet_t #m