% This code is for the Bonus Question

clear all;

% Model constants
radius = [0, 0.05, 0.085, 0.1]; % Include the center and all layer boundaries
r_tot = radius(end);  % in m
T_inf = 2; % oC
T_final = 5; % oC

% Material properties
k = [0.03, 0.6, 0.2];  
rho = [520, 1000, 1200];
cp = [1470, 4200, 3150];

% Spatial discretization
N = 41; 
dr = r_tot/(N-1);
r = [];

% Assign material properties to nodes
k_nodes = zeros(1, N);
rho_nodes = zeros(1, N);
cp_nodes = zeros(1, N);

layer_thickness = diff(radius);
total_thickness = sum(layer_thickness);
nodes_per_layer = round(N * layer_thickness / total_thickness);

% Apply material property to each node
index = 1; 
for layer = 1:length(nodes_per_layer)
    count = nodes_per_layer(layer);

    % Generate positions for nodes in this layer
    if layer == 1
        r_layer = linspace(radius(layer), radius(layer + 1), count); % First layer includes center
    else
        r_layer = linspace(radius(layer), radius(layer + 1), count + 1); % Include boundary so +1
        r_layer = r_layer(2:end); % Exclude duplicate boundary from previous layer
    end

    % Append node positions to r=[]
    r = [r, r_layer]; % Append radial positions

    % Assign material properties for this layer
    k_nodes(index:index + count - 1) = k(layer);
    rho_nodes(index:index + count - 1) = rho(layer);
    cp_nodes(index:index + count - 1) = cp(layer);

    % Update index for the next layer
    index = index + count;
end

r_len = length(r);

% Time discretization
dt = 5; % This is arbitrary because the calculated value did not converge
t_sim = 360000;
N_time = t_sim / dt + 1;

% Initial conditions
T_ini = 30; % oC

% Vary h from 50 to 180
h = 50:1:180; % Convection coefficients to test by the increment of 1 oC
time_to_cool = zeros(size(h)); % Array to store times for the graph

for h_index = 1:length(h)

    h_now = h(h_index); % Current convection coefficient

    T_old = ones(r_len, 1) * T_ini;
    T_new = T_old;
    Q_conv = 0;    % Reset convection energy
    A_outer = 4 * pi * r(end)^2; % Surface area of the cantaloupe

    % Time loop
    for i = 1:N_time - 1
        T_old = T_new; 
        
        % Center node
        T_new(1) = T_old(1) + (T_old(2) - T_old(1)) * (6 * k_nodes(1) * dt)/(rho_nodes(1) * cp_nodes(1) * dr^2);

        % Interior nodes
        for j = 2:N-1

            r_n = (j-1) * dr; % nodal position
            r_minus = r_n - dr / 2; 
            r_plus = r_n + dr / 2;

            T_new(j) = T_old(j) + ((T_old(j-1) - T_old(j)) * k_nodes(j-1)*r_minus^2 / dr + (T_old(j+1) - T_old(j)) * k_nodes(j+1) * r_plus^2 / dr ) / ((rho_nodes(j)*cp_nodes(j)*(r_n^3-r_minus^3) + rho_nodes(j)*cp_nodes(j)*(r_plus^3-r_n^3)) / (3*dt));

        end

        % Convection boundary node
        r_outer = r(end);
        r_inner = r(end) - dr / 2;

            T_new(N) = T_old(N) + (T_old(N-1) - T_old(N)) * 3 * k_nodes(N) * (r_outer - r_inner)^2 * dt / (dr * rho_nodes(N) * cp_nodes(N) * (r_outer^3 - (r_outer - dr/2)^3)) + (T_inf - T_old(N)) * 3 * h_now * r_outer^2 * dt / (rho_nodes(N) * cp_nodes(N) * (r_outer^3 - (r_outer - dr/2)^3));


        % Check stopping condition
        if T_new(nodes_per_layer(1)) <= T_final

            time_to_cool(h_index) = i * dt; % Record time for this `h` on the each node from 1 to 131

            break;
        end
    end
end

% Plot results
figure;
plot(h, time_to_cool / 3600, 'b-o', 'LineWidth', 2); % Convert time to hours
xlabel('Convection Coefficient, h [W/m^2-K]');
ylabel('Time to Reach 5°C at the first layer [hours]');
title('Bonus Question');
grid on;
