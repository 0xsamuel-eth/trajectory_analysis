clc; clear; close all;

% Parameters
rho = 1.225;    % Constant density (for now) (kg/m^3)
g = 9.81;       % Gravitational acceleration (m/s^2)
tb = 7;         % Burn time
mp = 42.5;      % Propellant mass (kg)
md = 192.4;     % Dead mass (kg)
mdot = mp/tb;   % Propellant mass flow rate (kg/s)
m = mp + md;    % Total mass before burn (kg)
Isp = 369;      % Specific impulse of the fuel (From what Arin found)
Thrust = Isp .* mdot .* g;       % Thrust force in N (constant for now)
Cd = 0.26;      % Drag coefficient
A = 0.25;       % Frontal area (m^2)
towerL = 5;     % Launch tower length (m)

%%

% Setup for solving
dt = 0.001;      % Time step in seconds
tBurn = 0:dt:tb;   % Time range for burn
tCoast = tb:dt:200;     % Time range for coast
t = [tBurn tCoast];     % Total time range

T = zeros(1, length(t));        % Initalize Thrust
D = zeros(1, length(t));        % Initalize Drag
Theta = zeros(1, length(t));    % Initalize Theta (Horizon Reference)
Theta(1) = 45;                  % Launch angle
M = zeros(1, length(t));        % Initalize Mass
M(1) = m;                       % Start mass is propellant plus deadweight
xdist = zeros(1, length(t));    % Intialize x distance traveled
ydist = zeros(1, length(t));    % Initialize y distance traveled
dist = zeros(1, length(t));     % Initialize distance traveled
x = zeros(1, length(t));        % Initalize x position
y = zeros(1, length(t));        % Initalize y position
Fn = zeros(1, length(t));       % Initialize normal force from the tower
Fx = zeros(1, length(t));       % Initalize force in x
Fy = zeros(1, length(t));       % Initalize force in y
Vx = zeros(1, length(t));       % Initialize velocity in x
Vy = zeros(1, length(t));       % Initialize velocity in y
Ax = zeros(1, length(t));       % Initialize acceleration in x
Ay = zeros(1, length(t));       % Initialize acceleration in y

% Rocket Equations for burn
for i = 2:length(tBurn)
    
    M(i) = m - mdot .* tBurn(i);
    
    if dist(i-1) <= towerL

        Fn(i) = M(i) * g * cosd(Theta(1));

    else

        Fn(i) = 0;

    end

    T(i) = Thrust * tBurn(i);  
    D(i) = 0.5 * Cd * rho * A * (Vx(i-1)^2 + Vy(i-1)^2); 
    Fx(i) = T(i) * cosd(Theta(i-1)) - D(i) * cosd(Theta(i-1)) ...
        - Fn(i) * sind(Theta(i-1));
    Fy(i) = T(i) * sind(Theta(i-1)) - D(i) * sind(Theta(i-1)) ...
        - (M(i) * g) + Fn(i) * cosd(Theta(i-1));
    Ax(i) = Fx(i) / M(i);
    Ay(i) = Fy(i) / M(i);
    Vx(i) = Vx(i-1) + Ax(i) * dt;
    Vy(i) = Vy(i-1) + Ay(i) * dt;
    Theta(i) = atand(Vy(i) / Vx(i));
    x(i) = x(i-1) + Vx(i) * dt;
    y(i) = y(i-1) + Vy(i) * dt;
    xdist(i) = xdist(i-1) + Vx(i) * dt;
    ydist(i) = ydist(i-1) + Vy(i) * dt;
    dist(i) = sqrt(xdist(i)^2 + ydist(i)^2);

end


%%

% Rocket Equations for coast
for i = length(tBurn) + 1:length(t)

    T(i) = 0; 
    D(i) = 0.5 * Cd * rho * A * (Vx(i-1)^2 + Vy(i-1)^2); 
    M(i) = md;
    Fn(i) = 0;
    Fx(i) = T(i) * cosd(Theta(i-1)) - D(i) * cosd(Theta(i-1)) - Fn(i) * sind(Theta(i-1));
    Fy(i) = T(i) * sind(Theta(i-1)) - D(i) * sind(Theta(i-1)) - (M(i) * g) + Fn(i) * cosd(Theta(i-1));
    Ax(i) = Fx(i) / M(i);
    Ay(i) = Fy(i) / M(i);
    Vx(i) = Vx(i-1) + Ax(i) * dt;
    Vy(i) = Vy(i-1) + Ay(i) * dt;
    Theta(i) = atand(Vy(i) / Vx(i));
    x(i) = x(i-1) + Vx(i) * dt;
    y(i) = y(i-1) + Vy(i) * dt;

end


V = sqrt(Vx.^2 + Vy.^2);
A = sqrt(Ax.^2 + Ay.^2);
F = sqrt(Fx.^2 + Fy.^2);



% Cut it off after the rocket hits the ground
V(y<0) = [];
A(y<0) = [];
F(y<0) = [];
T(y<0) = [];
D(y<0) = [];
M(y<0) = [];
Fx(y<0) = [];
Fy(y<0) = [];
Ax(y<0) = [];
Ay(y<0) = [];
Vx(y<0) = [];
Vy(y<0) = [];
x(y<0) = [];
Theta(y<0) = [];
t(y<0) = [];
y(y<0) = [];


%% Plotting

figure();
subplot (4,2,1);
plot(t, Fx);
xlabel('Time');
ylabel('Force (N)');
title('Force in X')

subplot (4,2,2);
plot(t, Fy);
xlabel('Time');
ylabel('Force (N)');
title('Force in Y')

subplot (4,2,3);
plot(t, Ax);
xlabel('Time');
ylabel('Acceleration (m/s^2)');
title('Acceleration in X')

subplot (4,2,4);
plot(t, Ay);
xlabel('Time');
ylabel('Acceleration (m/s^2)');
title('Acceleration in Y')

subplot (4,2,5);
plot(t, Vx);
xlabel('Time');
ylabel('Velocity (m/s)');
title('Velocity in X')

subplot (4,2,6);
plot(t, Vy);
xlabel('Time');
ylabel('Velocity (m/s)');
title('Velocity in Y')

subplot (4,2,7);
plot(t, x / 1000);
xlabel('Time');
ylabel('Position (km)');
title('Position in X')

subplot (4,2,8);
plot(t, y / 1000);
xlabel('Time');
ylabel('Position (km)');
title('Position in Y')

figure();
subplot(1,4,1);
plot(x / 1000, y / 1000);
xlabel('Range (km)');
ylabel('Height (km)');
title('2D Trajectory');

subplot(1,4,2);
plot(t, V);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Total Velocity');

subplot(1,4,3);
plot(t, Theta);
xlabel('Time (s)');
ylabel('Theta (deg)');
title('Horizon Reference Angle');

subplot(1,4,4);
plot(t, M);
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Mass');