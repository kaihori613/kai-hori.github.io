clear all;

% Model constants
radius = [0, 0.05, 0.085, 0.1]; % Include the center and all layer boundaries
r_tot = radius(end);  % in m
h = 50; % W/m2-K
T_inf = 2; % oC
T_final = 5; % oC

% Material properties
k = [0.03, 0.6, 0.2];  
rho = [520, 1000, 1200];
cp = [1470, 4200, 3150];

% Spatial discretization
N = 30; 
dr = r_tot/(N-1);
r = 0:dr:r_tot; 
r_len = length(r);

% Assign material properties to nodes
k_nodes = zeros(1, N);
rho_nodes = zeros(1, N);
cp_nodes = zeros(1, N);

layer_thickness = diff(radius);
total_thickness = sum(layer_thickness);
nodes_per_layer = round(N * layer_thickness / total_thickness);

% apply material property to each node
for layer = 1:length(radius) - 1
    index = find(r >= radius(layer) & r <= radius(layer + 1));
    k_nodes(index) = k(layer);
    rho_nodes(index) = rho(layer);
    cp_nodes(index) = cp(layer);
end

% Time discretization
%alpha_nodes = k_nodes ./ (rho_nodes .* cp_nodes);
%alpha_min = min(alpha_nodes(alpha_nodes > 0)); % Avoid zero values
%dt = 0.5 * dr^2 / alpha_min; % Stability criterion
dt = 20;
t_sim = 360000;
N_time = t_sim / dt + 1;
t = 0:dt:t_sim;

% Initial conditions
T_ini = 30; % oC
T_old = ones(r_len, 1) * T_ini;
T_new = T_old;
T_hist = zeros(r_len, N_time);

% Tolerance
tol = 10e-6; 

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
      
        T_new(j) = T_old(j) + ((T_old(j-1) - T_old(j)) * k_nodes(j)*r_minus^2 / dr + (T_old(j+1) - T_old(j)) * k_nodes(j) * r_plus^2 / dr ) / ((rho_nodes(j)*cp_nodes(j)*(r_n^3-r_minus^3) + rho_nodes(j)*cp_nodes(j)*(r_plus^3-r_n^3)) / (3*dt));
    end
    
    % Convection boundary node
    r_outer = r(end);
    r_inner = r(end) - dr / 2;

    T_new(N) = T_old(N) + (T_old(N-1) - T_old(N)) * 3 * k_nodes(N) * (r_outer - r_inner)^2 * dt / (dr * rho_nodes(N) * cp_nodes(N) * (r_outer^3 - (r_outer - dr/2)^3)) + (T_inf - T_old(N)) * 3 * h * r_outer^2 * dt / (rho_nodes(N) * cp_nodes(N) * (r_outer^3 - (r_outer - dr/2)^3));

    
    % Update history
    T_hist(:, i) = T_new;
   
    % Plot
    plot(r, T_new, 'b-', 'LineWidth', 2);
    xlabel('Position (r) [m]');
    ylabel('Temperature (T) [°C]');
    title(sprintf('Temperature Distribution at t = %.3f s', i * dt));
    grid on;
    pause(0.01);

    % Check stopping condition
    if T_new(15) < 5
        fprintf('Stopping simulation at time = %.2f seconds\n', i * dt);
        fprintf('Stopping simulation at time = %.2f hours\n', i * dt / 3600);
        fprintf('Temperature at r = 0.05m = %.3f oC\n', T_new(15));
        break;
    end

end

% Compute the volume for each node (control volume)
V = zeros(1, N); % Volume array
for i = 1:N
    if i == 1
        r_outside = r(i+1) - dr/2; % Half-shell for the first node
        V(i) = (4/3) * pi * r_outside^3;
    elseif i == N
        r_inside = r(i-1) + dr/2; % Half-shell for the last node
        V(i) = (4/3) * pi * (r(i)^3 - r_inside^3);
    else
        r_plus = r(i) + dr/2;
        r_minus = r(i) - dr/2;
        V(i) = (4/3) * pi * (r_plus^3 - r_minus^3);
    end
end

% Calculate the energy removed
Q_total = 0; % Initialize total energy
for i = 1:N
    d_T = T_ini - T_new(i); % Temperature difference
    Q_total = Q_total + rho_nodes(i) * cp_nodes(i) * V(i) * d_T; % Energy for each node
end

% Display the result
fprintf('Total energy required to cool the fruit: %.3f J\n', Q_total);