%% Description

% code to optimize cook time based on size, free-stream velocity, and air temperature of a 
% grill/oven for cooking 4X4X1.5 in steaks. given problem requires enegery
% use to be less than 20 kW and air temperature to be less than 200
% degrees celsius

clear; clc; close all;

%% defining constants

% steak (fat beef) properties
k_steak = 0.190;
alpha_steak = 0.13e-6;
rho_steak = 810;
cp_steak = 3540;
L_steak = 0.1016;
t_steak = 0.0381;

% temperatures
t_inf = 20;
t_desired = 66;

% grill
E_p = 0.25; % [kW/hr]
S_p = 1.7; % [$/kg]
grill_t = 0.01; % grill thickness
grill_l = L_steak + 0.01; % grill length
rho_steel = 7913;

%% Biot number coefficients table

Biot_table = [0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 ...
              2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 30.0 40.0 50.0 100.0 101; ...
              0.0998 0.1410 0.1987 0.2425 0.2791 0.3111 0.4328 0.5218 0.5932 ...
              0.6533 0.7051 0.7506 0.7910 0.8274 0.8603 1.0769 1.1925 1.2646 ...
              1.3138 1.3496 1.3766 1.3978 1.4149 1.4289 1.4961 1.5202 1.5325 ...
              1.5400 1.5552 1.5708; ...
              1.0017 1.0033 1.0066 1.0098 1.0130 1.0161 1.0311 1.0450 1.0580 ...
              1.0701 1.0814 1.0918 1.1016 1.1107 1.1191 1.1785 1.2102 1.2287 ...
              1.2403 1.2479 1.2532 1.2570 1.2598 1.2620 1.2699 1.2717 1.2723 ...
              1.2727 1.2731 1.2732];

%% read in air properties

air_prop = readmatrix("Air Properties.xlsx");

%% analysis using one term approximation of an infinite plate

% considered velocity and temperatures
V = linspace(0.01,5,100);
T = 100:200;

% preallocate array to store parameters
bigArr = zeros(length(V),length(T), 8);

for t_ind = 1:length(T)
    
    % calculate film temperature
    T_surf_avg = (20+150)/2;
    T_film = (T_surf_avg+T(t_ind))/2;
    
    % read in air properties based on film temp
    [~, temp_ind] = min(abs(air_prop(:,1) - T_film)); % index of closest film temp
    rho_air = air_prop(temp_ind,2);
    cp_air = air_prop(temp_ind,3);
    k_air = air_prop(temp_ind,4);
    nu_air = air_prop(temp_ind,5);
    Pr = air_prop(temp_ind,6);
    
    for v_ind = 1:length(V)

    Re = V(v_ind)*L_steak/nu_air; % calculate Reynolds number
    bigArr(v_ind, t_ind, 1) = Re; % store Re
    % calculate boundary layer thickness with 5mm of added buffer
    bt = 4.91*L_steak/Re^0.5 + 0.005;

    % calculate necessary gril size
    grill_w = 5*( 2*bt + L_steak ); % grill width
    grill_h = 2*(2*bt + t_steak); % grill height
    grill_A = grill_w*grill_h; % calculate grill frontal area
    
    % calculate steel required and cost
    S_area = 2*(grill_w*grill_h) + 2*(grill_w*grill_l) + 2*(grill_h*grill_l);
    S_w = S_area*grill_t*rho_steel;
    S_cost = S_w * S_p;
    bigArr(v_ind,t_ind,8) = S_cost;

    m_dot = rho_air * V(v_ind) * grill_A; % mass flow rate
    
    Q = m_dot*cp_air*(T(t_ind)-t_inf); % required energy input
    
if Q <= 20000 % start analysis if combination of V and T meets energy limit
    
    % calculate nusselts number
    if Re <= 50000 % if laminar
        Nu = 0.664 * Re^0.5 * Pr^(1/3);
    else
        Nu = 0.037 * Re^0.8 * Pr^(1/3); % turbulent
    end
    bigArr(v_ind,t_ind,3) = Nu; % store Nu
    
    h = Nu * k_air / L_steak; % calculate heat transfer coeff
    bigArr(v_ind, t_ind, 4) = h; % store h
    
    Bi = h * (t_steak/2) / k_steak; % biot number
    bigArr(v_ind, t_ind, 5) = Bi; % store Bi

    % read in coefficients based on closest biot number
    [~, bi_ind] = min(abs(Biot_table(1,:) - Bi)); % index of closest Biot number
    lam1 = Biot_table(2,bi_ind);
    A1 = Biot_table(3,bi_ind);
    
    if Bi > 0.1 % not lumped
        
        % calculate Fourier number
        tau = (1/-lam1)*log( ( (t_desired - T(t_ind))/(t_inf - T(t_ind) ) * (1/A1) ) );
        
        t_req = tau * (t_steak/2)^2 / alpha_steak; % solve for cook time
    
        t_surf = (20 - T(t_ind))*A1*exp(-lam1*tau)*cos(lam1) + T(t_ind); % check surface temp
        bigArr(v_ind, t_ind, 7) = t_surf; % store surface temp
    
        % nullify result if surface temperature exceeds 150
        if t_surf > 150
            t_req = NaN;
        bigArr(v_ind, t_ind, 7) = NaN; 
        end
        
    else % if lumped
        
        b = h/(rho_steak * cp_steak * (t_steak/2));
        t_req = (1/-b)*log( (t_desired - T(t_ind))/(t_inf - T(t_ind)) );
        bigArr(v_ind,t_ind,7) = t_req; % store surface temp
        
        % nullify result if surface temperature exceeds 150
        if t_req > 150
            t_req = NaN;
        end
        
    end
    
else
    % enter null result if Q exceeds energy limit
    Q = NaN;
    t_req = NaN;
end
    bigArr(v_ind, t_ind, 2) = Q; % store Q
    bigArr(v_ind, t_ind, 6) = t_req/60; % store cook time in minutes

    end
end

%% plot results

if length(V) > 1 && length(T) > 1 % if running entire sim
    
% setup and plot surface plot
figure
[vel,temp] = meshgrid(T,V);
surf(vel,temp,bigArr(:,:,6),bigArr(:,:,2), 'EdgeColor','none')
% plot formatting
fs = 12; col = colorbar;
ylabel(col,'Energy Required [W]','FontSize',fs,'FontWeight','bold', 'Rotation',270);
col.Label.Position(1) = 3;
xlabel("Air Temperature [C]",'FontSize',fs,'FontWeight','bold')
ylabel("Airflow Velocity [m/s]",'FontSize',fs,'FontWeight','bold')
zlabel("Cook Time [Minutes]",'FontSize',fs,'FontWeight','bold')
view([5, 5, 5])

figure
surf(vel,temp,bigArr(:,:,6),bigArr(:,:,8), 'EdgeColor','none')
% plot formatting
fs = 12; col = colorbar;
ylabel(col,'Steel Cost [$]','FontSize',fs,'FontWeight','bold', 'Rotation',270);
col.Label.Position(1) = 3;
xlabel("Air Temperature [C]",'FontSize',fs,'FontWeight','bold')
ylabel("Airflow Velocity [m/s]",'FontSize',fs,'FontWeight','bold')
zlabel("Cook Time [Minutes]",'FontSize',fs,'FontWeight','bold')
view([5, 5, 5])

end

%% best results

min_t = min(bigArr(:,:,6),[],'all');
[v_op, t_op] = find(bigArr(:,:,6)==min_t);
Re_op = bigArr(v_op(1),t_op(1),1);
Q_op = bigArr(v_op(1),t_op(1),2);
h_op = bigArr(v_op(1),t_op(1),4);
Ts_op = bigArr(v_op(1),t_op(1),7);
S_c = bigArr(v_op(1),t_op(1),8);
E_cost = Q_op/1000*E_p;





