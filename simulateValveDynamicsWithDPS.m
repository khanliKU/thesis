function result = simulateValveDynamicsWithDPS(...
var_min_air_gap,...
var_c_r        ,...
var_c_a        ,...
var_h_C        ,...
var_r_Cout     ,...
var_D_w        ,...
var_t_cfr      ,...
var_t_cfa      ,...
var_t_mba      ,...
var_t_mta      ,...
var_r_A        ,...
var_r_cham     ,...
var_h_cham     ,...
var_t_shell    ,...
var_PR         ,...
var_CSF        ,...
var_h_air      ,...
var_mu         ,...
var_mu0        ,...
var_CWCC       ,...
var_erc        ,...
dt             ,...
tao_multiplier ,...
P_U            ,...
T_U            ,...
P_D            ,...
T_D            ,...
DPS_vol        ,...
DPS_orifice_A  ,...
writeFile      ,...
fileName       ,...
showPlot        ...
)
load 'analytics.mat'

load_b_h_430f;

Nturns               = floor(func_N(var_D_w,var_PR,var_c_r,var_h_C,var_r_A,var_r_Cout,var_t_cfr))
func_coil_length(var_D_w, var_PR, var_c_r, var_h_C, var_r_A, var_r_Cout, var_t_cfr)
resistance = func_R(var_D_w,var_PR,var_c_r,var_h_C,var_r_A,var_r_Cout,var_erc,var_t_cfr)
max_I = func_I_max(var_CSF,var_CWCC,var_D_w)
Nturns*max_I
Vdd = resistance * max_I
max_H = Nturns * max_I / var_h_C


% check if it event opens
MR_total = func_R_total(var_c_a,var_c_r,var_h_C,var_h_air,var_h_cham,var_mu,var_mu0,var_r_A,var_r_Cout,var_r_cham,var_t_cfa,var_t_mba,var_t_mta,var_t_shell);
flux = Nturns * max_I / MR_total;
mag_field = Nturns * max_I / var_h_C;
mag_flux_density_cham  = flux /(pi * (var_r_A^2 - var_r_cham^2)/cos(atan(var_h_cham/var_r_cham)));
mag_flux_density_arm  = flux /(pi * (var_r_A - var_r_cham)^2);
mag_flux_density_arm_cham  = flux /(pi * var_r_A^2);
mag_flux_density_shell  = flux /(pi*((var_c_r+var_r_Cout+var_t_shell)^2-(var_c_r+var_r_Cout)^2));
mag_flux_density_limit  = correctFlux(mag_field ,bh430F);
if max([mag_flux_density_arm_cham , mag_flux_density_shell ]) > mag_flux_density_limit 
    mag_flux_density  = mag_flux_density_limit ;
else
    mag_flux_density  = max([mag_flux_density_arm_cham , mag_flux_density_shell ]);
end
F_mag_max  = 0.5 * (mag_flux_density ^2 * var_r_A^2 * pi) / var_mu0
%

tao_open   = tao_multiplier * (Nturns^2/func_R_total(var_c_a,var_c_r,var_h_C,var_h_air,var_h_cham,var_mu,var_mu0,var_r_A,var_r_Cout,var_r_cham,var_t_cfa,var_t_mba,var_t_mta,var_t_shell))/resistance
tao_close  = tao_multiplier  * (Nturns^2/func_R_total(var_c_a,var_c_r,var_h_C,var_min_air_gap,var_h_cham,var_mu,var_mu0,var_r_A,var_r_Cout,var_r_cham,var_t_cfa,var_t_mba,var_t_mta,var_t_shell))/resistance
t_max = ceil((tao_open+tao_close)/dt) % s
signal(1,1:t_max) = 1;
t(1,1:t_max) = (0:(t_max-1))*dt; % s

voltage(1,1:t_max) = Vdd;
voltage_L(1,1:t_max) = voltage(1);
current(1,1:t_max) = 0;
MR_total(1,1:t_max) = func_R_total(var_c_a,var_c_r,var_h_C,var_h_air,var_h_cham,var_mu,var_mu0,var_r_A,var_r_Cout,var_r_cham,var_t_cfa,var_t_mba,var_t_mta,var_t_shell);
flux(1,1:t_max) = 0;
inductance(1,1:t_max) = Nturns^2 / MR_total(1);
di_dt(1,1:t_max) = voltage_L(1) / inductance(1);
mag_flux_density_cham(1,1:t_max) = 0;
mag_flux_density_arm(1,1:t_max) = 0;
mag_flux_density_arm_cham(1,1:t_max) = 0;
mag_flux_density_shell(1,1:t_max) = 0;
mag_flux_density_limit(1,1:t_max) = 0;
mag_flux_density(1,1:t_max) = 0;
mag_field(1,1:t_max) = 0;


pilot(1,1:t_max) = 0; % 1 = open, 0 = closed

gam = 1.4;
R_gc = 296.8; % J/kgK
Cv = 0.743; % kJ/kg.K
polythropic_index = gam ;

P_u(1,1:t_max) = P_U; % bar
P_d(1,1:t_max) = P_D; % bar
P_cham(1,1:t_max) = P_u(1); % bar
P_airgap(1,1:t_max) = P_u(1); % bar
P_dps(1,1:t_max) = P_d(1); % bar

T_u(1,1:t_max) = T_U; % K
T_d(1,1:t_max) = T_D; % K
T_cham(1,1:t_max) = T_u(1); % K
T_airgap(1,1:t_max) = T_u(1); % K
T_dps(1,1:t_max) = T_U; % K

m_dot_u_cham(1,1:t_max) = 0;
m_dot_cham_dps(1,1:t_max) = 0;
m_dot_dps_d(1,1:t_max) = 0;
m_dot_cham_airgap(1,1:t_max) = 0;
m_dot_cham_net(1,1:t_max) = 0;

x0 = 0; %m
x_pre = 0.5E-3; %m
x(1,1:t_max) = x0; % m
v(1,1:t_max) = 0; % m/s
a(1,1:t_max) = 0; % m/s^2

Ks = 25000; % N/m
density430F = 7800;
W_armature = double(subs(W_magArmature,[h_C,r_A,rho_permeable,t_cfa],[var_h_C,var_r_A,density430F,var_t_cfa]))
%{
air gap and chamber volumes are depending on the armature displacement.
%}
r_outlet = 2 * var_h_air;
r_inlet =  r_outlet * sqrt(2);

func_V_airgap = matlabFunction(h_air * (r_A + c_r)^2 * pi + 0.25 * pi * h_C * ((r_A + c_r)^2 - r_A^2) );
func_V_chamber = matlabFunction(5E-3 * (r_A + c_r)^2 * pi - h_air * r_A^2 * pi + 0.25 * pi * h_C * ((r_A + c_r)^2 - r_A^2));
Vol_airgap(1,1:t_max) = func_V_airgap(var_c_r,var_h_C,var_h_air,var_r_A);
Vol_chamber(1,1:t_max) = func_V_chamber(var_c_r,var_h_C,var_h_air,var_r_A);

M_airgap(1,1:t_max) = P_airgap(1) * Vol_airgap(1) / (T_airgap(1) * R_gc);
M_cham(1,1:t_max) = P_cham(1) * Vol_chamber(1) / (T_cham(1) * R_gc);
M_dps(1,1:t_max) = P_dps(1) * DPS_vol / (T_dps(1) * R_gc);


Cd_inout = 0.6;
Cd_clearance = 0.5;

A_inlet = pi * r_inlet^2; % * Cd_inout;
A_outlet = pi * r_outlet^2; % * Cd_inout;

A_clearance = ((var_r_A + var_c_r)^2 - var_r_A^2 )* pi; % * Cd_clearance;

F_airgap(1,1:t_max) = P_airgap(1) * var_r_A^2 * pi;
F_cham(1,1:t_max) = P_cham(1) * (var_r_A^2-r_outlet^2) * pi;
F_spring(1,1:t_max) = (x_pre + x(1)) * Ks;
F_mag(1,1:t_max) = 0;
F(1,1:t_max) = - F_airgap(1) + F_cham(1) - F_spring(1);
force_to_overcome = F(1)
if (F_mag_max * 0.8) < abs(force_to_overcome)
    return
end

f = waitbar(0, 'Starting');
last_iteration = 1;
max_v = 0;
min_v = 0;
close_valve = 0;
%tic
%try

for i = 2:t_max

    if mod(i,round(t_max/100)) == 0
        waitbar(i/t_max, f, sprintf('Progress: %.1f %% - %.1f seconds elapsed.\n%.1f seconds projected - %.1f seconds remain.', i/t_max*100,toc,toc/(i/t_max),toc/(i/t_max)-toc));
    end
    
    % Motion with last steps forces

    x(i) = checkImag(x(i-1) + v(i-1)*dt + 0.5*a(i-1)*dt^2); % m
    v(i) = checkImag(v(i-1) + a(i-1)*dt); % m/s
    a(i) = checkImag(F(i-1) / W_armature); % m/s^2
    if x(i) <= 0 && a(i) <= 0
        x(i) = 0; % m
        v(i) = 0; % m/s
        a(i) = 0; % m/s^2
    elseif (var_h_air - x(i)) <= var_min_air_gap && a(i) >= 0
        x(i) = var_h_air - var_min_air_gap; % m
        v(i) = 0; % m/s
        a(i) = 0; % m/s^2
    end

    signal(i) = signal(i-1);

    % Old Design Magnetism
    current(i) = current(i-1) + di_dt(i-1) * dt;

    % flux(i) = func_phi(var_D_w,current(i),0.95,)
    MR_total(i) = func_R_total(var_c_a,var_c_r,var_h_C,var_h_air - x(i),var_h_cham,var_mu,var_mu0,var_r_A,var_r_Cout,var_r_cham,var_t_cfa,var_t_mba,var_t_mta,var_t_shell);
    flux(i) = Nturns * current(i) / MR_total(i);
    mag_field(i) = Nturns * current(i) / var_h_C;
    mag_flux_density_cham(i) = flux(i)/(pi * (var_r_A^2 - var_r_cham^2)/cos(atan(var_h_cham/var_r_cham)));
    mag_flux_density_arm(i) = flux(i)/(pi * (var_r_A - var_r_cham)^2);
    mag_flux_density_arm_cham(i) = flux(i)/(pi * var_r_A^2);
    mag_flux_density_shell(i) = flux(i)/(pi*((var_c_r+var_r_Cout+var_t_shell)^2-(var_c_r+var_r_Cout)^2));
    mag_flux_density_limit(i) = correctFlux(mag_field(i),bh430F);
    if max([mag_flux_density_arm_cham(i), mag_flux_density_shell(i)]) > mag_flux_density_limit(i)
        mag_flux_density(i) = mag_flux_density_limit(i);
    else
        mag_flux_density(i) = max([mag_flux_density_arm_cham(i), mag_flux_density_shell(i)]);
    end
    if current(i) ~= 0
        inductance(i) = Nturns * (mag_flux_density(i) * var_r_A^2 * pi) / current(i);
    else
        inductance(i) = Nturns^2 / MR_total(i);
    end
    F_mag(i) = 0.5 * (mag_flux_density(i)^2 * var_r_A^2 * pi) / var_mu0;
    if (signal(i) > 0)
        applied_V = Vdd;
    else
        applied_V = 0;
    end
    voltage_L(i) = applied_V - current(i) * resistance;
    di_dt(i) = voltage_L(i) / inductance(i) ;

    % Fluids
    Vol_airgap(i) = func_V_airgap(var_c_r,var_h_C,var_h_air - x(i),var_r_A);
    Vol_chamber(i) = func_V_chamber(var_c_r,var_h_C,var_h_air - x(i),var_r_A);
    
        % Adiabatic Compression
    [T_cham(i), P_cham(i)] = volChange(Vol_chamber(i-1),Vol_chamber(i),T_cham(i-1),P_cham(i-1),polythropic_index);

    [T_airgap(i), P_airgap(i)] = volChange(Vol_airgap(i-1),Vol_airgap(i),T_airgap(i-1),P_airgap(i-1),1);

        % Outlet Flows
    [P_dps(i),T_dps(i),M_dps(i),P_d(i),T_d(i),dummy,m_dot_dps_d(i)] = ...
        gasExchange(P_dps(i-1),T_cham(i-1),M_dps(i-1),DPS_vol,P_D,T_D,inf,inf,gam,DPS_orifice_A,R_gc,Cv,dt,m_dot_dps_d(i-1));
    
    [P_cham(i),T_cham(i),M_cham(i),P_dps(i),T_dps(i),M_dps(i),m_dot_cham_dps(i)] = ...
        gasExchange(P_cham(i),T_cham(i),M_cham(i-1),Vol_chamber(i),P_dps(i),T_dps(i),M_dps(i),DPS_vol,gam,A_outlet*(x(i)/var_h_air),R_gc,Cv,dt,m_dot_cham_dps(i-1));

%    if P_cham(i)<P_airgap(i)
%        ahgjhgjh =1;
%    end
%gasExchange(P1,T1,M1,vol1,P2,T2,M2,vol2,gam,A,R,Cv,dt,mdot_old)
        % Inlet Flows
    [P_u(i),T_u(i),dummy,P_cham(i),T_cham(i),M_cham(i),m_dot_u_cham(i)] = ...
        gasExchange(P_U,T_U,inf,inf,P_cham(i),T_cham(i),M_cham(i),Vol_chamber(i),gam,A_inlet,R_gc,Cv,dt,m_dot_u_cham(i-1));

        % Internal Flow
    [P_cham(i),T_cham(i),M_cham(i),P_airgap(i),T_airgap(i),M_airgap(i),m_dot_cham_airgap(i)] = ...
        gasExchange(P_cham(i),T_cham(i),M_cham(i),Vol_chamber(i),P_airgap(i),T_airgap(i),M_airgap(i-1),Vol_airgap(i),gam,A_clearance,R_gc,Cv,dt,m_dot_cham_airgap(i-1));
    
    %M_airgap(i) = M_airgap(i-1);
    %m_dot_cham_airgap(i) = m_dot_cham_airgap(i-1);

    m_dot_cham_net(i) = m_dot_u_cham(i) - m_dot_cham_dps(i) - m_dot_cham_airgap(i);
        %{
    [P_temp, T_temp, M_temp] = tank_discharge_io(...
        P_cham(i),T_cham(i),M_cham(i-1),...
        [T_U,     T_airgap(i), T_D],...
        [P_U,     P_airgap(i), P_D],...
        [A_inlet, A_clearance, A_outlet*(x(i)/var_h_air)],...
        gam,R_gc,dt); % bar, K, kg
    P_cham(i) = checkImag(P_temp);
    T_cham(i) = checkImag(T_temp);
    M_cham(i) = checkImag(M_temp);

    [P_temp, T_temp, M_temp] = tank_discharge_io(...
        P_airgap(i),T_airgap(i),M_airgap(i-1),...
        [T_cham(i)],...
        [P_cham(i)],...
        [A_clearance],...
        gam,R_gc,dt); % bar, K, kg

    P_airgap(i) = checkImag(P_temp);
    T_airgap(i) = checkImag(T_temp);
    M_airgap(i) = checkImag(M_temp);
        %}

    % Forces

    F_airgap(i) = P_airgap(i) * var_r_A^2 * pi;
    F_cham(i) = P_cham(i) * (var_r_A^2 - r_outlet^2) * pi;
    F_spring(i) = (x_pre + x(i)) * Ks;
    F(i) = checkImag(- F_airgap(i) + F_cham(i) - F_spring(i) + F_mag(i));
    
    last_iteration = i;
end
%{
catch
    errordlg('Imaginary number!','Computation Error');
    i
end
%}
close(f);
%toc
max_Force = max(F)
orifice_opening_t = t(find(x == (var_h_air - var_min_air_gap),1))
if isempty(orifice_opening_t)
    orifice_opening_t = 0
    orifice_closing_t = 0
    total_flow = 0
else
    tao_open = orifice_opening_t
    orifice_closing_t = t(find(v == min(v),1)) - tao_open
    total_flow  = sum(m_dot_cham_dps)*dt
end
total_power = sum(current*Vdd)*dt
weight      = func_weight(var_PR, var_c_r, var_h_C, var_r_A, var_r_Cout, var_t_cfa, var_t_cfr, var_t_mba, var_t_mta, var_t_shell)

save_struct.var_min_air_gap      = var_min_air_gap      ;
save_struct.var_c_r              = var_c_r              ;
save_struct.var_c_a              = var_c_a              ;
save_struct.var_h_C              = var_h_C              ;
save_struct.var_r_Cout           = var_r_Cout           ;
save_struct.var_D_w              = var_D_w              ;
save_struct.var_t_cfr            = var_t_cfr            ;
save_struct.var_t_cfa            = var_t_cfa            ;
save_struct.var_t_mba            = var_t_mba            ;
save_struct.var_t_mta            = var_t_mta            ;
save_struct.var_r_A              = var_r_A              ;
save_struct.var_r_cham           = var_r_cham           ;
save_struct.var_h_cham           = var_h_cham           ;
save_struct.var_t_shell          = var_t_shell          ;
save_struct.var_PR               = var_PR               ;
save_struct.var_CSF              = var_CSF              ;
save_struct.var_h_air            = var_h_air            ;
save_struct.var_mu               = var_mu               ;
save_struct.var_mu0              = var_mu0              ;
save_struct.var_CWCC             = var_CWCC             ;
save_struct.var_erc              = var_erc              ;
save_struct.Nturns               = Nturns               ;
save_struct.P                    = P_U                  ;
save_struct.Ks                   = Ks                   ;
save_struct.orifice_opening_t    = orifice_opening_t    ;
save_struct.orifice_closing_t    = orifice_closing_t    ;
save_struct.total_flow           = total_flow           ;
save_struct.total_power          = total_power          ;
save_struct.weight               = weight               ;
save_struct.deltaP               = P_U - P_dps(end)     ;
save_struct.mdot                 = m_dot_cham_dps(end)  ;
save_struct.chamber_rho          = M_cham(end)/Vol_chamber(end);
save_struct.nozzle_A             = DPS_orifice_A;
          
if writeFile
    writetable(struct2table(save_struct),fileName,"WriteMode",'append')
end

result.save_struct = save_struct;

% Reduce data size
if length(t) == t_max
    %tic
    reduce_factor = 100;
    result.pilot = downsample(movmean(pilot(1:last_iteration),reduce_factor),reduce_factor);
    result.signal = downsample(movmean(signal(1:last_iteration),reduce_factor),reduce_factor);
    result.t = downsample(t(1:last_iteration),reduce_factor);

    result.P_airgap = Pa_to_B(downsample(movmean(P_airgap(1:last_iteration),reduce_factor),reduce_factor));
    result.P_cham = Pa_to_B(downsample(movmean(P_cham(1:last_iteration),reduce_factor),reduce_factor));
    result.P_u = Pa_to_B(downsample(movmean(P_u(1:last_iteration),reduce_factor),reduce_factor));
    result.P_d = Pa_to_B(downsample(movmean(P_d(1:last_iteration),reduce_factor),reduce_factor));
    result.P_dps = Pa_to_B(downsample(movmean(P_dps(1:last_iteration),reduce_factor),reduce_factor));


    result.T_airgap = downsample(movmean(T_airgap(1:last_iteration),reduce_factor),reduce_factor);
    result.T_cham = downsample(movmean(T_cham(1:last_iteration),reduce_factor),reduce_factor);
    result.T_u = downsample(movmean(T_u(1:last_iteration),reduce_factor),reduce_factor);
    result.T_d = downsample(movmean(T_d(1:last_iteration),reduce_factor),reduce_factor);
    result.T_dps = downsample(movmean(T_dps(1:last_iteration),reduce_factor),reduce_factor);

    result.M_airgap = downsample(movmean(M_airgap(1:last_iteration),reduce_factor),reduce_factor);
    result.M_cham = downsample(movmean(M_cham(1:last_iteration),reduce_factor),reduce_factor);
    result.M_dps = downsample(movmean(M_dps(1:last_iteration),reduce_factor),reduce_factor);

    result.Vol_airgap = downsample(movmean(Vol_airgap(1:last_iteration),reduce_factor),reduce_factor);
    result.Vol_chamber = downsample(movmean(Vol_chamber(1:last_iteration),reduce_factor),reduce_factor);

    result.m_dot_dps_d = downsample(movmean(m_dot_dps_d(1:last_iteration),reduce_factor),reduce_factor);    
    result.m_dot_cham_dps = downsample(movmean(m_dot_cham_dps(1:last_iteration),reduce_factor),reduce_factor);                               
    result.m_dot_cham_airgap = downsample(movmean(m_dot_cham_airgap(1:last_iteration),reduce_factor),reduce_factor);
    result.m_dot_u_cham = downsample(movmean(m_dot_u_cham(1:last_iteration),reduce_factor),reduce_factor);
    result.m_dot_cham_net = downsample(movmean(m_dot_cham_net(1:last_iteration),reduce_factor),reduce_factor);

    result.di_dt = downsample(movmean(di_dt(1:last_iteration),reduce_factor),reduce_factor);
    result.current = downsample(movmean(current(1:last_iteration),reduce_factor),reduce_factor);
    result.voltage_L = downsample(movmean(voltage_L(1:last_iteration),reduce_factor),reduce_factor);
    result.voltage = downsample(movmean(voltage(1:last_iteration),reduce_factor),reduce_factor);
    result.inductance = downsample(movmean(inductance(1:last_iteration),reduce_factor),reduce_factor);
    result.MR_total = downsample(movmean(MR_total(1:last_iteration),reduce_factor),reduce_factor);
    result.flux = downsample(movmean(flux(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_field = downsample(movmean(mag_field(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_flux_density = downsample(movmean(mag_flux_density(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_flux_density_cham = downsample(movmean(mag_flux_density_cham(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_flux_density_arm = downsample(movmean(mag_flux_density_arm(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_flux_density_shell = downsample(movmean(mag_flux_density_shell(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_flux_density_limit = downsample(movmean(mag_flux_density_limit(1:last_iteration),reduce_factor),reduce_factor);
    result.mag_flux_density_arm_cham = downsample(movmean(mag_flux_density_arm_cham(1:last_iteration),reduce_factor),reduce_factor);

    result.x = downsample(movmean(x(1:last_iteration),reduce_factor),reduce_factor);
    result.v = downsample(movmean(v(1:last_iteration),reduce_factor),reduce_factor);
    result.a = downsample(movmean(a(1:last_iteration),reduce_factor),reduce_factor);
    
    result.F = downsample(movmean(F(1:last_iteration),reduce_factor),reduce_factor);
    result.F_mag = downsample(movmean(F_mag(1:last_iteration),reduce_factor),reduce_factor);
    result.F_spring = downsample(movmean(F_spring(1:last_iteration),reduce_factor),reduce_factor);
    result.F_cham = downsample(movmean(F_cham(1:last_iteration),reduce_factor),reduce_factor);
    result.F_airgap = downsample(movmean(F_airgap(1:last_iteration),reduce_factor),reduce_factor);
    %toc
end

if showPlot

    pilot = downsample(movmean(pilot(1:last_iteration),reduce_factor),reduce_factor);
    signal = downsample(movmean(signal(1:last_iteration),reduce_factor),reduce_factor);
    t = downsample(t(1:last_iteration),reduce_factor);
    P_airgap = Pa_to_B(downsample(movmean(P_airgap(1:last_iteration),reduce_factor),reduce_factor));
    P_cham = Pa_to_B(downsample(movmean(P_cham(1:last_iteration),reduce_factor),reduce_factor));
    P_u = Pa_to_B(downsample(movmean(P_u(1:last_iteration),reduce_factor),reduce_factor));
    P_d = Pa_to_B(downsample(movmean(P_d(1:last_iteration),reduce_factor),reduce_factor));
    P_dps = Pa_to_B(downsample(movmean(P_dps(1:last_iteration),reduce_factor),reduce_factor));
    T_airgap = downsample(movmean(T_airgap(1:last_iteration),reduce_factor),reduce_factor);
    T_cham = downsample(movmean(T_cham(1:last_iteration),reduce_factor),reduce_factor);
    T_u = downsample(movmean(T_u(1:last_iteration),reduce_factor),reduce_factor);
    T_d = downsample(movmean(T_d(1:last_iteration),reduce_factor),reduce_factor);
    T_dps = downsample(movmean(T_dps(1:last_iteration),reduce_factor),reduce_factor);
    M_airgap = downsample(movmean(M_airgap(1:last_iteration),reduce_factor),reduce_factor);
    M_cham = downsample(movmean(M_cham(1:last_iteration),reduce_factor),reduce_factor);
    M_dps = downsample(movmean(M_dps(1:last_iteration),reduce_factor),reduce_factor);
    Vol_airgap = downsample(movmean(Vol_airgap(1:last_iteration),reduce_factor),reduce_factor);
    Vol_chamber = downsample(movmean(Vol_chamber(1:last_iteration),reduce_factor),reduce_factor);
    m_dot_dps_d = downsample(movmean(m_dot_dps_d(1:last_iteration),reduce_factor),reduce_factor);    
    m_dot_cham_dps = downsample(movmean(m_dot_cham_dps(1:last_iteration),reduce_factor),reduce_factor);                               
    m_dot_cham_airgap = downsample(movmean(m_dot_cham_airgap(1:last_iteration),reduce_factor),reduce_factor);
    m_dot_u_cham = downsample(movmean(m_dot_u_cham(1:last_iteration),reduce_factor),reduce_factor);
    m_dot_cham_net = downsample(movmean(m_dot_cham_net(1:last_iteration),reduce_factor),reduce_factor);
    di_dt = downsample(movmean(di_dt(1:last_iteration),reduce_factor),reduce_factor);
    current = downsample(movmean(current(1:last_iteration),reduce_factor),reduce_factor);
    voltage_L = downsample(movmean(voltage_L(1:last_iteration),reduce_factor),reduce_factor);
    voltage = downsample(movmean(voltage(1:last_iteration),reduce_factor),reduce_factor);
    inductance = downsample(movmean(inductance(1:last_iteration),reduce_factor),reduce_factor);
    MR_total = downsample(movmean(MR_total(1:last_iteration),reduce_factor),reduce_factor);
    flux = downsample(movmean(flux(1:last_iteration),reduce_factor),reduce_factor);
    mag_field = downsample(movmean(mag_field(1:last_iteration),reduce_factor),reduce_factor);
    mag_flux_density = downsample(movmean(mag_flux_density(1:last_iteration),reduce_factor),reduce_factor);
    mag_flux_density_cham = downsample(movmean(mag_flux_density_cham(1:last_iteration),reduce_factor),reduce_factor);
    mag_flux_density_arm = downsample(movmean(mag_flux_density_arm(1:last_iteration),reduce_factor),reduce_factor);
    mag_flux_density_shell = downsample(movmean(mag_flux_density_shell(1:last_iteration),reduce_factor),reduce_factor);
    mag_flux_density_limit = downsample(movmean(mag_flux_density_limit(1:last_iteration),reduce_factor),reduce_factor);
    mag_flux_density_arm_cham = downsample(movmean(mag_flux_density_arm_cham(1:last_iteration),reduce_factor),reduce_factor);
    x = downsample(movmean(x(1:last_iteration),reduce_factor),reduce_factor);
    v = downsample(movmean(v(1:last_iteration),reduce_factor),reduce_factor);
    a = downsample(movmean(a(1:last_iteration),reduce_factor),reduce_factor);
    F = downsample(movmean(F(1:last_iteration),reduce_factor),reduce_factor);
    F_mag = downsample(movmean(F_mag(1:last_iteration),reduce_factor),reduce_factor);
    F_spring = downsample(movmean(F_spring(1:last_iteration),reduce_factor),reduce_factor);
    F_cham = downsample(movmean(F_cham(1:last_iteration),reduce_factor),reduce_factor);
    F_airgap = downsample(movmean(F_airgap(1:last_iteration),reduce_factor),reduce_factor);

orifice_opening_t = t(find(v == max(v),1))
tao_open = orifice_opening_t;
orifice_closing_t = t(find(v == min(v),1)) - tao_open
%tao_open = tao_open *1e3;
figure()
plot(t,m_dot_dps_d,t,m_dot_cham_dps,t,m_dot_cham_airgap,t,m_dot_u_cham),xline(tao_open)
legend('DPS -> D','cham -> DPS','cham -> airgap','U -> cham')

figure()
plot(t,P_airgap,t,P_cham,t,P_u,t,P_d,t,P_dps),xline(tao_open),legend('Airgap','Chamber','U','D','DPS','Power off')
title('Pressures'),grid on,xlabel('Time [ms]'),ylabel('[Bar]')
figure()
subplot(2,1,1)
plot(t,P_airgap,t,P_cham,t,P_dps),xline(tao_open),legend('Airgap','Chamber','DPS','Power off')
title('Pressures'),grid on,xlabel('Time [ms]'),ylabel('[Bar]')
subplot(2,1,2)
plot(t,T_airgap,t,T_cham,t,T_dps),xline(tao_open),legend('Airgap','Chamber','DPS','Power off')
title('Temperatures'),grid on,xlabel('Time [ms]'),ylabel('[K]')
figure() 
subplot(4,2,1)
plot(t,current),xline(tao_open)
title('Current'),grid on,xlabel('Time [ms]'),ylabel('[A]')
subplot(4,2,2)
plot(t,voltage_L),xline(tao_open)
title('Induction Voltage'),grid on,xlabel('Time [ms]'),ylabel('[V]')
subplot(4,2,3)
plot(t,MR_total),xline(tao_open)
title('Magnetic Reluctance'),grid on,xlabel('Time [ms]'),ylabel('[H^{-1}]'),xlim([0 inf])
subplot(4,2,4)
plot(t,inductance),xline(tao_open)
title('Inductance'),grid on,xlabel('Time [ms]'),ylabel('[H]')
subplot(4,2,5)
plot(t,flux),xline(tao_open)
title('Magnetic Flux'),grid on,xlabel('Time [ms]'),ylabel('[Wb]')
subplot(4,2,6)
plot(t,current.*Vdd),xline(tao_open)
title('Electrical Power'),grid on,xlabel('Time [ms]'),ylabel('[W]')
subplot(4,2,7)
plot(t,mag_flux_density),xline(tao_open)
title('Magnetic Flux Density'),grid on,xlabel('Time [ms]'),ylabel('[T]')
subplot(4,2,8)
plot(t,mag_field),xline(tao_open)
title('Magnetic Field'),grid on,xlabel('Time [ms]'),ylabel('[H]')
load_b_h_430f;
figure()
plot(mag_field*1e-3,mag_flux_density,bh430F.H_Oe*1E3/(4*pi)*1e-3,bh430F.B_Gauss/1E4),legend('B','Material Limit')
title('B - H Curve'),grid on,xlabel('Magnetic Field Strength H [ kA/m ]'),ylabel('Magnetic Flux Density B [ T ]')
figure()
semilogy(t,M_airgap,t,M_cham,t,M_dps),xline(tao_open),legend('Airgap','Chamber','DPS','Power off')
title('Masses'),grid on,xlabel('Time [ms]'),ylabel('[kg]')
figure()
subplot(2,1,1)
semilogy(t,Vol_airgap,t,Vol_chamber),xline(tao_open),legend('Airgap','Chamber','Power off')
title('Volumes'),grid on,xlabel('Time [ms]'),ylabel('[m^3]')
subplot(2,1,2)
semilogy(t,M_airgap./Vol_airgap,t,M_cham./Vol_chamber),xline(tao_open),legend('Airgap','Chamber','Power off')
title('Density'), grid on,xlabel('Time [ms]'),ylabel('[kg/m^3]')
figure()
subplot(4,1,1)
plot(t,m_dot_cham_net),title('Net'),grid on,xlabel('Time [ms]'),ylabel('[kg/s]'),xline(tao_open)
subplot(4,1,2)
plot(t,m_dot_cham_dps),title('Chamber -> DPS'),grid on,xlabel('Time [ms]'),ylabel('[kg/s]'),xline(tao_open)
subplot(4,1,3)
plot(t,m_dot_cham_airgap),title('Chamber -> Airgap'),grid on,xlabel('Time [ms]'),ylabel('[kg/s]'),xline(tao_open)
subplot(4,1,4)
plot(t,m_dot_u_cham),title('U -> Chamber'),grid on,xlabel('Time [ms]'),ylabel('[kg/s]'),xline(tao_open)
figure()
subplot(3,1,1)
plot(t,x*1e3)
title('Translation'),grid on,xlabel('Time [ms]'),ylabel('[mm]'),xline(tao_open)
subplot(3,1,2)
plot(t,v)
title('Velocity'),grid on,xlabel('Time [ms]'),ylabel('[m/s]'),xline(tao_open)
subplot(3,1,3)
plot(t,a)
title('Acceleration'),grid on,xlabel('Time [ms]'),ylabel('[m/s^2]'),xline(tao_open)
figure()
plot(t,F,t,F_mag,t,-F_spring,t,F_cham,t,-F_airgap),xline(tao_open),legend('Sum','Mag','Spring','Cham','Airgap','Power off')
title('Forces'),grid on,xlabel('Time [ms]'),ylabel('[N]')
figure()
plot(t,F,t,F_mag,t,-F_spring),xline(tao_open),legend('Sum','Mag','Spring','Power off')
title('Forces'),grid on,xlabel('Time [ms]'),ylabel('[N]')
figure()
semilogy(t,mag_flux_density_limit,t,mag_flux_density_shell,t,mag_flux_density_arm,t,mag_flux_density_cham,t,mag_flux_density_arm_cham,t,mag_flux_density),xline(tao_open)
legend('Limit','Shell','Armature','Chamfer','Arm-Cham','Used','Power off')
title('Magnetic Flux Densities'),grid on,xlabel('Time [ms]'),ylabel('[T]')
figure()
subplot(4,1,1)
plot(t,current),xline(tao_open)
title('Current'),grid on,xlabel('Time [ms]'),ylabel('[A]')
subplot(4,1,2)
plot(t,F,t,F_mag,t,-F_spring),xline(tao_open),legend('Sum','Mag','Spring','Power off')
title('Forces'),grid on,xlabel('Time [ms]'),ylabel('[N]')
subplot(4,1,3)
plot(t,x*1e3)
title('Translation'),grid on,xlabel('Time [ms]'),ylabel('[mm]'),xline(tao_open)
subplot(4,1,4)
plot(t,m_dot_cham_dps),title('Chamber -> DPS'),grid on,xlabel('Time [ms]'),ylabel('[kg/s]'),xline(tao_open)
end

end