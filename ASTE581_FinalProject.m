%% Initialize Problem and Constants
clear;

global mu;
mu_Earth = 398600; % km^3/s^2 (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
mu = mu_Earth;
global Re;
Re = 6378.137; % km (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
global J2;
J2 = 1082.63e-6; % Earth (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
T_GEO = 23.9345*3600; % sec (1 sidereal day) (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
global i_GEO;
i_GEO = 0; % rad
global r_GEO;
r_GEO = (mu*(T_GEO/(2*pi))^2)^(1/3); % km
global v_GEO;
v_GEO = sqrt(mu/r_GEO); % km/s
a_SE = 149.598e6; % km (Sun-Earth) (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
global mu_s;
mu_s = 132712e6; % km^3/s^2 (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
v_SE = sqrt(mu_s/a_SE); %km/s
f0_SE = 0; % rad
i_SE = 0; % rad (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
a_EM = 0.3844e6; %km (Earth-Moon) (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html)
global mu_m;
mu_m = 0.00490e6; % km^3/s^2 (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html)
v_EM = sqrt(mu_Earth/a_EM); % km/s
f0_EM = 0; % rad
i_EM = 5.145 * pi/180; % rad (referenced https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html)
n_svs = 34; % Number of SVs in constellation
fstep = 2*pi/n_svs; % rad (True anomaly increments)

% Initial state
X0 = zeros(18*n_svs, 1);
for i = 1:n_svs
    
    % Earth-Satellite
    f = (i - 1)*fstep;
    r0_esat = [r_GEO*cos(f)*cos(i_GEO) r_GEO*sin(f)*cos(i_GEO) r_GEO*sin(i_GEO)]'; % km
    v0_esat = [-v_GEO*sin(f)*cos(i_GEO) v_GEO*cos(f)*cos(i_GEO) v_GEO*sin(i_GEO)]'; % km/s
    X0((18*(i - 1) + 1):(18*(i - 1) + 3), 1) = r0_esat;
    X0((18*(i - 1) + 4):(18*(i - 1) + 6), 1) = v0_esat;
    
    % Sun-Earth
    r0_se = [a_SE*cos(f0_SE)*cos(i_SE) a_SE*sin(f0_SE)*cos(i_SE) a_SE*sin(i_SE)]'; % km
    v0_se = [-v_SE*sin(f0_SE)*cos(i_SE) v_SE*cos(f0_SE)*cos(i_SE) v_SE*sin(i_SE)]'; % km/s
    X0((18*(i - 1) + 7):(18*(i - 1) + 9), 1) = r0_se;
    X0((18*(i - 1) + 10):(18*(i - 1) + 12), 1) = v0_se;
    
    % Earth-Moon
    r0_em = [a_EM*cos(f0_EM)*cos(i_EM) a_EM*sin(f0_EM)*cos(i_EM) a_EM*sin(i_EM)]'; % km
    v0_em = [-v_EM*sin(f0_EM)*cos(i_EM) v_EM*cos(f0_EM)*cos(i_EM) v_EM*sin(i_EM)]'; % km/s
    X0((18*(i - 1) + 13):(18*(i - 1) + 15), 1) = r0_se;
    X0((18*(i - 1) + 16):(18*(i - 1) + 18), 1) = v0_se;
    
end

% Time span
tspan = 0:T_GEO:365*T_GEO; % sec
tspan_plot = 0:T_GEO/20:365*T_GEO; % sec

% Maneuver threshold
d_RMP = 1200; % km (diameter of RMP Earth coverage spot)
circum_Earth = 2*pi*Re; % km (circumfrence of Earth)
theta_threshold = 2*pi/(circum_Earth/d_RMP); % rad (angular separation)
SVa = [r_GEO*cos(0) r_GEO*sin(0) 0];
SVb = [r_GEO*cos(theta_threshold) r_GEO*sin(theta_threshold) 0];
threshold = norm(SVb - SVa); % km (separation distance threshold)

% Function handle
fn = @ode_J2SM_cowells;
global effects;

% Tolerances for ODE solver
options = odeset('RelTol', 1.0e-13, 'InitialStep', 1.0e-12, 'AbsTol', 1.0e-10);  

% Plot constants
orbitcolors = {'r', 'g', 'b', 'c', 'm'}; % Colors for orbit plane plots
t = 0:2*pi/360:2*pi; % Circular cross-section of Earth
e1 = zeros(1, length(t));
e2 = zeros(1, length(t));
for j = 1:length(t)
    e1(j) = cos(t(j));
    e2(j) = sin(t(j));
end
r = linspace(0, 1, n_svs); % Colors for time-series plots
b = [linspace(0.5, 1, n_svs/2) flip(linspace(0.5, 1, n_svs/2))];
g = flip(linspace(0, 1, n_svs));
colors = [r', b', g'];
legend_labels = cell(1, n_svs); % Labels for time-series legends
for i = 1:n_svs
    if i < 10
        legend_labels{i} = ['SV0' num2str(i)];
    else
        legend_labels{i} = ['SV' num2str(i)];
    end
end

%% Two-Body Problem
tic;
disp('Two-Body Problem');

% Integrate 3D two-body problem
effects = [0 0 0];
X_twobody = zeros(length(tspan), 18*n_svs);
for i = 1:n_svs
    [t, X_twobody(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers
[X_twobody, maneuvercount_twobody, maneuveridx_twobody] ...
    = GEO_maneuver_detect(X_twobody, tspan, n_svs, threshold, fstep, fn, options);

% Integrate 3D two-body problem for plane plot
X_twobody_plot = zeros(length(tspan_plot), 18*n_svs);
for i = 1:n_svs
    [t, X_twobody_plot(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan_plot, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers for plane plot
[X_twobody_plot, maneuvercount_twobody_plot, maneuveridx_twobody_plot] ...
    = GEO_maneuver_detect(X_twobody_plot, tspan_plot, n_svs, threshold, fstep, fn, options);

% Plot x-y plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobody_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobody_plot(maneuveridx_twobody_plot(i):maneuveridx_twobody_plot(i+1), 1)/Re, ...
        X_twobody_plot(maneuveridx_twobody_plot(i):maneuveridx_twobody_plot(i+1), 2)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Y (R_{Earth})');
title({'Two Body Problem', 'X-Y Plane', ['Maneuvers Conducted: ' num2str(maneuvercount_twobody_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBody_XYPlane.png');
close(gcf);

% Plot x-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobody_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobody_plot(maneuveridx_twobody_plot(i):maneuveridx_twobody_plot(i+1), 1)/Re, ...
        X_twobody_plot(maneuveridx_twobody_plot(i):maneuveridx_twobody_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem', 'X-Z Plane', ['Maneuvers Conducted: ' num2str(maneuvercount_twobody_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBody_XZPlane.png');
close(gcf);

% Plot y-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobody_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobody_plot(maneuveridx_twobody_plot(i):maneuveridx_twobody_plot(i+1), 2)/Re, ...
        X_twobody_plot(maneuveridx_twobody_plot(i):maneuveridx_twobody_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('Y (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem', 'Y-Z Plane', ['Maneuvers Conducted: ' num2str(maneuvercount_twobody_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBody_YZPlane.png');
close(gcf);

% Plot separation distance
figure;
distOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        svA = X_twobody(i, (18*(j - 1) + 1):(18*(j - 1) + 3));
        if j < n_svs
            svB = X_twobody(i, (18*(j) + 1):(18*(j) + 3));
        else
            svB = X_twobody(i, 1:3);
        end
        distOverTime(i, j) = norm(svA - svB);
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, distOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
plot(tspan/T_GEO, ones(1, length(distOverTime(:, 1)))*threshold, 'r--');
ax = gca;
lh = ax.Legend;
lh.String = lh.String(1:n_svs);
hold off;
xlabel('Time (Orbital Period)');
ylabel('Distance (km)');
title({'Two Body Problem', 'Separation Distance to N+1 SV', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobody)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBody_SepDist.png');
close(gcf);

% Plot orbital radius
figure;
rOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        rOverTime(i, j) = norm(X_twobody(i, (18*(j - 1) + 1):(18*(j - 1) + 3)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, rOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time (Orbital Period)');
ylabel('Orbital Radius (km)');
title({'Two Body Problem', 'Orbital Radius', ['Maneuvers Conducted: ' num2str(maneuvercount_twobody)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBody_Position.png');
close(gcf);

% Plot orbital velocity
figure;
vOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        vOverTime(i, j) = norm(X_twobody(i, (18*(j - 1) + 4):(18*(j - 1) + 6)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, vOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
xlabel('Time (Orbital Period)');
ylabel('Orbital Velocity (km/s)');
title({'Two Body Problem', 'Orbital Velocity', ['Maneuvers Conducted: ' num2str(maneuvercount_twobody)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBody_Velocity.png');
close(gcf);

toc;
%% Two-Body Problem (J2)
tic;
disp('Two-Body Problem (J2)');

% Integrate 3D two-body problem with J2 perturbations
effects = [1 0 0];
X_twobodyJ2 = zeros(length(tspan), 18*n_svs);
for i = 1:n_svs
    [t, X_twobodyJ2(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers
[X_twobodyJ2, maneuvercount_twobodyJ2, maneuveridx_twobodyJ2] ...
    = GEO_maneuver_detect(X_twobodyJ2, tspan, n_svs, threshold, fstep, fn, options);

% Integrate 3D two-body problem with J2 perturbations for plane plot
X_twobodyJ2_plot = zeros(length(tspan_plot), 18*n_svs);
for i = 1:n_svs
    [t, X_twobodyJ2_plot(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan_plot, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers for plane plot
[X_twobodyJ2_plot, maneuvercount_twobodyJ2_plot, maneuveridx_twobodyJ2_plot] ...
    = GEO_maneuver_detect(X_twobodyJ2_plot, tspan_plot, n_svs, threshold, fstep, fn, options);

% Plot x-y plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2_plot(maneuveridx_twobodyJ2_plot(i):maneuveridx_twobodyJ2_plot(i+1), 1)/Re, ...
        X_twobodyJ2_plot(maneuveridx_twobodyJ2_plot(i):maneuveridx_twobodyJ2_plot(i+1), 2)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Y (R_{Earth})');
title({'Two Body Problem (J2)', 'X-Y Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_XYPlane.png');
close(gcf);

% Plot x-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2_plot(maneuveridx_twobodyJ2_plot(i):maneuveridx_twobodyJ2_plot(i+1), 1)/Re, ...
        X_twobodyJ2_plot(maneuveridx_twobodyJ2_plot(i):maneuveridx_twobodyJ2_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem (J2)', 'X-Z Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_XZPlane.png');
close(gcf);

% Plot y-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2_plot(maneuveridx_twobodyJ2_plot(i):maneuveridx_twobodyJ2_plot(i+1), 2)/Re, ...
        X_twobodyJ2_plot(maneuveridx_twobodyJ2_plot(i):maneuveridx_twobodyJ2_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('Y (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem (J2)', 'Y-Z Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_YZPlane.png');
close(gcf);

% Plot separation distance
figure;
distOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        svA = X_twobodyJ2(i, (18*(j - 1) + 1):(18*(j - 1) + 3));
        if j < n_svs
            svB = X_twobodyJ2(i, (18*(j) + 1):(18*(j) + 3));
        else
            svB = X_twobodyJ2(i, 1:3);
        end
        distOverTime(i, j) = norm(svA - svB);
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, distOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
plot(tspan/T_GEO, ones(1, length(distOverTime(:, 1)))*threshold, 'r--');
ax = gca;
lh = ax.Legend;
lh.String = lh.String(1:n_svs);
hold off;
xlabel('Time (Orbital Period)');
ylabel('Distance (km)');
title({'Two Body Problem (J2)', 'Separation Distance to N+1 SV', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_SepDist.png');
close(gcf);

% Plot orbital radius
figure;
rOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        rOverTime(i, j) = norm(X_twobodyJ2(i, (18*(j - 1) + 1):(18*(j - 1) + 3)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, rOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time (Orbital Period)');
ylabel('Orbital Radius (km)');
title({'Two Body Problem (J2)', 'Orbital Radius', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_Position.png');
close(gcf);

% Plot orbital velocity
figure;
vOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        vOverTime(i, j) = norm(X_twobodyJ2(i, (18*(j - 1) + 4):(18*(j - 1) + 6)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, vOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
xlabel('Time (Orbital Period)');
ylabel('Orbital Velocity (km/s)');
title({'Two Body Problem (J2)', 'Orbital Velocity', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_Velocity.png');
close(gcf);

% Plot difference between two body (J2) and two body problems - position
figure;
xStateDiff = zeros(length(tspan), n_svs);
yStateDiff = zeros(length(tspan), n_svs);
zStateDiff = zeros(length(tspan), n_svs);
for i = 1:n_svs
    xStateDiff(:, i) = X_twobodyJ2(:, (18*(i - 1) + 1)) - X_twobody(:, (18*(i - 1) + 1));
    yStateDiff(:, i) = X_twobodyJ2(:, (18*(i - 1) + 2)) - X_twobody(:, (18*(i - 1) + 2));
    zStateDiff(:, i) = X_twobodyJ2(:, (18*(i - 1) + 3)) - X_twobody(:, (18*(i - 1) + 3));
end
plot(tspan/T_GEO, mean(xStateDiff, 2));
hold on;
plot(tspan/T_GEO, mean(yStateDiff, 2));
plot(tspan/T_GEO, mean(zStateDiff, 2));
hold off;
legend({'x', 'y', 'z'}, 'location', 'northwest');
xlabel('Time (Orbital Period)');
ylabel('Mean Difference in Position (km)');
title({'Two Body Problem (J2) and Two Body Problem', ...
    'Position Components Comparison', ...
    ['Difference in Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2 - maneuvercount_twobody)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_diffTwoBodyPos.png');
close(gcf);

% Plot difference between two body (J2) and two body problems - velocity
figure;
vxStateDiff = zeros(length(tspan), n_svs);
vyStateDiff = zeros(length(tspan), n_svs);
vzStateDiff = zeros(length(tspan), n_svs);
for i = 1:n_svs
    vxStateDiff(:, i) = X_twobodyJ2(:, (18*(i - 1) + 4)) - X_twobody(:, (18*(i - 1) + 4));
    vyStateDiff(:, i) = X_twobodyJ2(:, (18*(i - 1) + 5)) - X_twobody(:, (18*(i - 1) + 5));
    vzStateDiff(:, i) = X_twobodyJ2(:, (18*(i - 1) + 6)) - X_twobody(:, (18*(i - 1) + 6));
end
plot(tspan/T_GEO, mean(vxStateDiff, 2));
hold on;
plot(tspan/T_GEO, mean(vyStateDiff, 2));
plot(tspan/T_GEO, mean(vzStateDiff, 2));
hold off;
legend({'vx', 'vy', 'vz'}, 'location', 'northwest');
xlabel('Time (Orbital Period)');
ylabel('Mean Difference in Velocity (km/s)');
title({'Two Body Problem (J2) and Two Body Problem', ...
    'Velocity Components Comparison', ...
    ['Difference in Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2 - maneuvercount_twobody)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2_diffTwoBodyVel.png');
close(gcf);

toc;
%% Two-Body Problem (J2, Sun-Earth)
tic;
disp('Two-Body Problem (J2, S-E)');

% Integrate 3D two-body problem with J2 and Sun-Earth perturbations
effects = [1 1 0];
X_twobodyJ2S = zeros(length(tspan), 18*n_svs);
for i = 1:n_svs
    [t, X_twobodyJ2S(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers
[X_twobodyJ2S, maneuvercount_twobodyJ2S, maneuveridx_twobodyJ2S] ...
    = GEO_maneuver_detect(X_twobodyJ2S, tspan, n_svs, threshold, fstep, fn, options);

% Integrate 3D two-body problem with J2 and Sun-Earth perturbations for
% plane plot
X_twobodyJ2S_plot = zeros(length(tspan_plot), 18*n_svs);
for i = 1:n_svs
    [t, X_twobodyJ2S_plot(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan_plot, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers for plane plot
[X_twobodyJ2S_plot, maneuvercount_twobodyJ2S_plot, maneuveridx_twobodyJ2S_plot] ...
    = GEO_maneuver_detect(X_twobodyJ2S_plot, tspan_plot, n_svs, threshold, fstep, fn, options);

% Plot x-y plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2S_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2S_plot(maneuveridx_twobodyJ2S_plot(i):maneuveridx_twobodyJ2S_plot(i+1), 1)/Re, ...
        X_twobodyJ2S_plot(maneuveridx_twobodyJ2S_plot(i):maneuveridx_twobodyJ2S_plot(i+1), 2)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Y (R_{Earth})');
title({'Two Body Problem (J2, S-E)', 'X-Y Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_XYPlane.png');
close(gcf);

% Plot x-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2S_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2S_plot(maneuveridx_twobodyJ2S_plot(i):maneuveridx_twobodyJ2S_plot(i+1), 1)/Re, ...
        X_twobodyJ2S_plot(maneuveridx_twobodyJ2S_plot(i):maneuveridx_twobodyJ2S_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem (J2, S-E)', 'X-Z Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_XZPlane.png');
close(gcf);

% Plot y-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2S_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2S_plot(maneuveridx_twobodyJ2S_plot(i):maneuveridx_twobodyJ2S_plot(i+1), 2)/Re, ...
        X_twobodyJ2S_plot(maneuveridx_twobodyJ2S_plot(i):maneuveridx_twobodyJ2S_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('Y (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem (J2, S-E)', 'Y-Z Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_YZPlane.png');
close(gcf);

% Plot separation distance
figure;
distOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        svA = X_twobodyJ2S(i, (18*(j - 1) + 1):(18*(j - 1) + 3));
        if j < n_svs
            svB = X_twobodyJ2S(i, (18*(j) + 1):(18*(j) + 3));
        else
            svB = X_twobodyJ2S(i, 1:3);
        end
        distOverTime(i, j) = norm(svA - svB);
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, distOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
plot(tspan/T_GEO, ones(1, length(distOverTime(:, 1)))*threshold, 'r--');
ax = gca;
lh = ax.Legend;
lh.String = lh.String(1:n_svs);
hold off;
xlabel('Time (Orbital Period)');
ylabel('Distance (km)');
title({'Two Body Problem (J2, S-E)', 'Separation Distance to N+1 SV', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_SepDist.png');
close(gcf);

% Plot orbital radius
figure;
rOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        rOverTime(i, j) = norm(X_twobodyJ2S(i, (18*(j - 1) + 1):(18*(j - 1) + 3)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, rOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time (Orbital Period)');
ylabel('Orbital Radius (km)');
title({'Two Body Problem (J2, S-E)', 'Orbital Radius', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_Position.png');
close(gcf);

% Plot orbital velocity
figure;
vOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        vOverTime(i, j) = norm(X_twobodyJ2S(i, (18*(j - 1) + 4):(18*(j - 1) + 6)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, vOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
xlabel('Time (Orbital Period)');
ylabel('Orbital Velocity (km/s)');
title({'Two Body Problem (J2, S-E)', 'Orbital Velocity', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_Velocity.png');
close(gcf);

% Plot difference between two body (J2, S-E) and two body (J2) problems - position
figure;
xStateDiff = zeros(length(tspan), n_svs);
yStateDiff = zeros(length(tspan), n_svs);
zStateDiff = zeros(length(tspan), n_svs);
for i = 1:n_svs
    xStateDiff(:, i) = X_twobodyJ2S(:, (18*(i - 1) + 1)) - X_twobodyJ2(:, (18*(i - 1) + 1));
    yStateDiff(:, i) = X_twobodyJ2S(:, (18*(i - 1) + 2)) - X_twobodyJ2(:, (18*(i - 1) + 2));
    zStateDiff(:, i) = X_twobodyJ2S(:, (18*(i - 1) + 3)) - X_twobodyJ2(:, (18*(i - 1) + 3));
end
plot(tspan/T_GEO, mean(xStateDiff, 2));
hold on;
plot(tspan/T_GEO, mean(yStateDiff, 2));
plot(tspan/T_GEO, mean(zStateDiff, 2));
hold off;
legend({'x', 'y', 'z'}, 'location', 'northwest');
xlabel('Time (Orbital Period)');
ylabel('Mean Difference in Position (km)');
title({'Two Body Problem (J2, S-E) and Two Body Problem (J2)', ...
    'Position Components Comparison', ...
    ['Difference in Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S - maneuvercount_twobodyJ2)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_diffTwoBodyJ2Pos.png');
close(gcf);

% Plot difference between two body (J2, S-E) and two body (J2) problems - velocity
figure;
vxStateDiff = zeros(length(tspan), n_svs);
vyStateDiff = zeros(length(tspan), n_svs);
vzStateDiff = zeros(length(tspan), n_svs);
for i = 1:n_svs
    vxStateDiff(:, i) = X_twobodyJ2S(:, (18*(i - 1) + 4)) - X_twobodyJ2(:, (18*(i - 1) + 4));
    vyStateDiff(:, i) = X_twobodyJ2S(:, (18*(i - 1) + 5)) - X_twobodyJ2(:, (18*(i - 1) + 5));
    vzStateDiff(:, i) = X_twobodyJ2S(:, (18*(i - 1) + 6)) - X_twobodyJ2(:, (18*(i - 1) + 6));
end
plot(tspan/T_GEO, mean(vxStateDiff, 2));
hold on;
plot(tspan/T_GEO, mean(vyStateDiff, 2));
plot(tspan/T_GEO, mean(vzStateDiff, 2));
hold off;
legend({'vx', 'vy', 'vz'}, 'location', 'northwest');
xlabel('Time (Orbital Period)');
ylabel('Mean Difference in Velocity (km/s)');
title({'Two Body Problem (J2, S-E) and Two Body Problem (J2)', ...
    'Velocity Components Comparison', ...
    ['Difference in Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2S - maneuvercount_twobodyJ2)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2S_diffTwoBodyJ2Vel.png');
close(gcf);

toc;
%% Two Body Problem (J2, Sun-Earth, Earth-Moon)
tic;
disp('Two Body Problem (J2, S-E, E-M)');

% Integrate 3D two-body problem with J2, Sun-Earth, and Earth-Moon perturbations
effects = [1 1 1];
X_twobodyJ2SM = zeros(length(tspan), 18*n_svs);
for i = 1:n_svs
    [t, X_twobodyJ2SM(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers
[X_twobodyJ2SM, maneuvercount_twobodyJ2SM, maneuveridx_twobodyJ2SM] ...
    = GEO_maneuver_detect(X_twobodyJ2SM, tspan, n_svs, threshold, fstep, fn, options);

% Integrate 3D two-body problem with J2, Sun-Earth, and Earth-Moon
% perturbations for plane plot
X_twobodyJ2SM_plot = zeros(length(tspan_plot), 18*n_svs);
for i = 1:n_svs
    [t, X_twobodyJ2SM_plot(:, (18*(i - 1) + 1):(18*(i - 1) + 18))] ...
        = ode45(fn, tspan_plot, X0((18*(i - 1) + 1):(18*(i - 1) + 18), 1), options);
end

% Detect maneuvers for plane plot
[X_twobodyJ2SM_plot, maneuvercount_twobodyJ2SM_plot, maneuveridx_twobodyJ2SM_plot] ...
    = GEO_maneuver_detect(X_twobodyJ2SM_plot, tspan_plot, n_svs, threshold, fstep, fn, options);

% Plot x-y plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2SM_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2SM_plot(maneuveridx_twobodyJ2SM_plot(i):maneuveridx_twobodyJ2SM_plot(i+1), 1)/Re, ...
        X_twobodyJ2SM_plot(maneuveridx_twobodyJ2SM_plot(i):maneuveridx_twobodyJ2SM_plot(i+1), 2)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Y (R_{Earth})');
title({'Two Body Problem (J2, S-E, E-M)', 'X-Y Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_XYPlane.png');
close(gcf);

% Plot x-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2SM_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2SM_plot(maneuveridx_twobodyJ2SM_plot(i):maneuveridx_twobodyJ2SM_plot(i+1), 1)/Re, ...
        X_twobodyJ2SM_plot(maneuveridx_twobodyJ2SM_plot(i):maneuveridx_twobodyJ2SM_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('X (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem (J2, S-E, E-M)', 'X-Z Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_XZPlane.png');
close(gcf);

% Plot y-z plane
figure;
plot(e1, e2, 'k');
hold on;
for i = 1:length(maneuveridx_twobodyJ2SM_plot)-1 % Plot each orbit between maneuvers
    coloridx = mod(i, length(orbitcolors));
    if coloridx == 0
        coloridx = length(orbitcolors);
    end
    plot(X_twobodyJ2SM_plot(maneuveridx_twobodyJ2SM_plot(i):maneuveridx_twobodyJ2SM_plot(i+1), 2)/Re, ...
        X_twobodyJ2SM_plot(maneuveridx_twobodyJ2SM_plot(i):maneuveridx_twobodyJ2SM_plot(i+1), 3)/Re, ...
        orbitcolors{coloridx});
    hold on;
end
grid on;
xlabel('Y (R_{Earth})');
ylabel('Z (R_{Earth})');
title({'Two Body Problem (J2, S-E, E-M)', 'Y-Z Plane', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM_plot)]});
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_YZPlane.png');
close(gcf);

% Plot separation distance
figure;
distOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        svA = X_twobodyJ2SM(i, (18*(j - 1) + 1):(18*(j - 1) + 3));
        if j < n_svs
            svB = X_twobodyJ2SM(i, (18*(j) + 1):(18*(j) + 3));
        else
            svB = X_twobodyJ2SM(i, 1:3);
        end
        distOverTime(i, j) = norm(svA - svB);
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, distOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
plot(tspan/T_GEO, ones(1, length(distOverTime(:, 1)))*threshold, 'r--');
ax = gca;
lh = ax.Legend;
lh.String = lh.String(1:n_svs);
hold off;
xlabel('Time (Orbital Period)');
ylabel('Distance (km)');
title({'Two Body Problem (J2, S-E, E-M)', 'Separation Distance to N+1 SV', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_SepDist.png');
close(gcf);

% Plot orbital radius
figure;
rOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        rOverTime(i, j) = norm(X_twobodyJ2SM(i, (18*(j - 1) + 1):(18*(j - 1) + 3)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, rOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('Time (Orbital Period)');
ylabel('Orbital Radius (km)');
title({'Two Body Problem (J2, S-E, E-M)', 'Orbital Radius', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_Position.png');
close(gcf);

% Plot orbital velocity
figure;
vOverTime = zeros(length(tspan), n_svs);
for i = 1:length(tspan)
    for j = 1:n_svs
        vOverTime(i, j) = norm(X_twobodyJ2SM(i, (18*(j - 1) + 4):(18*(j - 1) + 6)));
    end
end
for i = 1:n_svs
    plot(tspan/T_GEO, vOverTime(:, i), ...
        'Color', colors(i, :));
    hold on;
end
hold off;
legend(legend_labels, 'Location', 'northeastoutside', 'NumColumns', 2);
xlabel('Time (Orbital Period)');
ylabel('Orbital Velocity (km/s)');
title({'Two Body Problem (J2, S-E, E-M)', 'Orbital Velocity', ...
    ['Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_Velocity.png');
close(gcf);

% Plot difference between two body (J2, S-E, E-M) and two body (J2, S-E) problems - position
figure;
xStateDiff = zeros(length(tspan), n_svs);
yStateDiff = zeros(length(tspan), n_svs);
zStateDiff = zeros(length(tspan), n_svs);
for i = 1:n_svs
    xStateDiff(:, i) = X_twobodyJ2SM(:, (18*(i - 1) + 1)) - X_twobodyJ2S(:, (18*(i - 1) + 1));
    yStateDiff(:, i) = X_twobodyJ2SM(:, (18*(i - 1) + 2)) - X_twobodyJ2S(:, (18*(i - 1) + 2));
    zStateDiff(:, i) = X_twobodyJ2SM(:, (18*(i - 1) + 3)) - X_twobodyJ2S(:, (18*(i - 1) + 3));
end
plot(tspan/T_GEO, mean(xStateDiff, 2));
hold on;
plot(tspan/T_GEO, mean(yStateDiff, 2));
plot(tspan/T_GEO, mean(zStateDiff, 2));
hold off;
legend({'x', 'y', 'z'}, 'location', 'northwest');
xlabel('Time (Orbital Period)');
ylabel('Mean Difference in Position (km)');
title({'Two Body Problem (J2, S-E, E-M) and Two Body Problem (J2, S-E)', ...
    'Position Components Comparison', ...
    ['Difference in Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM - maneuvercount_twobodyJ2S)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_diffTwoBodyJ2SPos.png');
close(gcf);

% Plot difference between two body (J2, S-E, E-M) and two body (J2, S-E) problems - velocity
figure;
vxStateDiff = zeros(length(tspan), n_svs);
vyStateDiff = zeros(length(tspan), n_svs);
vzStateDiff = zeros(length(tspan), n_svs);
for i = 1:n_svs
    vxStateDiff(:, i) = X_twobodyJ2SM(:, (18*(i - 1) + 4)) - X_twobodyJ2S(:, (18*(i - 1) + 4));
    vyStateDiff(:, i) = X_twobodyJ2SM(:, (18*(i - 1) + 5)) - X_twobodyJ2S(:, (18*(i - 1) + 5));
    vzStateDiff(:, i) = X_twobodyJ2SM(:, (18*(i - 1) + 6)) - X_twobodyJ2S(:, (18*(i - 1) + 6));
end
plot(tspan/T_GEO, mean(vxStateDiff, 2));
hold on;
plot(tspan/T_GEO, mean(vyStateDiff, 2));
plot(tspan/T_GEO, mean(vzStateDiff, 2));
hold off;
legend({'vx', 'vy', 'vz'}, 'location', 'northwest');
xlabel('Time (Orbital Period)');
ylabel('Mean Difference in Velocity (km/s)');
title({'Two Body Problem (J2, S-E, E-M) and Two Body Problem (J2, S-E)', ...
    'Velocity Components Comparison', ...
    ['Difference in Maneuvers Conducted: ' num2str(maneuvercount_twobodyJ2SM - maneuvercount_twobodyJ2S)]});
grid on;
saveas(gcf, 'ASTE581\FinalProject\Figures\ASTE581_FinalProject_TwoBodyJ2SM_diffTwoBodyJ2SVel.png');
close(gcf);

toc;
%% Validation - DirecTV-15
tic;
disp('Validation - DirecTV-15');

% Initial parameters (4/18/2022 16:00:00) (NASA JPL HORIZONS)
r0_se = [-1.322574708254664E+08 -7.120535511008461E+07 4.289022599551827E+03]; % Sun-Earth
v0_se = [1.362658131963260E+01 -2.633474893703986E+01 2.439047330110355E-03];
r0_em = [-2.177334626275463E+05 -2.938494594337748E+05 -6.810797786048497E+02]; % Earth-Moon
v0_em = [8.699532350602149E-01 -6.268144526653440E-01 -9.832388893488264E-02];

% Time span
tspan = 0:60:seconds(datetime(2022, 4, 28, 16, 0, 0) - datetime(2022, 4, 18, 16, 0, 0)); % sec

% Integrate 3D two-body problem
disp('Two Body Problem');
% Initial parameters (4/18/2022 16:00:00) (STK)
r0_esat = [40504.609381 -11709.621832 10.658693]; % Earth-Satellite
v0_esat = [0.853874 2.953818 -0.000352];
X0 = [r0_esat'; v0_esat'; r0_se'; v0_se'; r0_em'; v0_em'];
% Final parameters (4/28/2022 16:00:00) (STK)
rf_esat = [41902.455145 -4681.303450 9.687367];
vf_esat = [0.341356 3.055759 -0.000479];
% Propagate orbit
effects = [0 0 0];
[~, validate_twobody] = ode45(fn, tspan, X0, options);
% Print results
diffpos_twobody = norm(validate_twobody(end, 1:3) - rf_esat);
disp(['Difference in Position = ' num2str(diffpos_twobody) ' km']);
diffvel_twobody = norm(validate_twobody(end, 4:6) - vf_esat);
disp(['Difference in Velovity = ' num2str(diffvel_twobody) ' km/s']);

% Integrate 3D two-body problem with J2 perturbations
disp('Two Body Problem (J2)');
% Initial parameters (4/18/2022 16:00:00) (STK)
r0_esat = [40509.080078 -11694.145720 10.655925]; % Earth-Satellite
v0_esat = [0.852746 2.954144 -0.000352];
X0 = [r0_esat'; v0_esat'; r0_se'; v0_se'; r0_em'; v0_em'];
% Final parameters (4/28/2022 16:00:00) (STK)
rf_esat = [41925.624709 -4469.032710 9.637240];
vf_esat = [0.325876 3.057448 -0.000484];
% Propagate orbit
effects = [1 0 0];
[~, validate_twobodyJ2] = ode45(fn, tspan, X0, options);
% Print results
diffpos_twobodyJ2 = norm(validate_twobodyJ2(end, 1:3) - rf_esat);
disp(['Difference in Position = ' num2str(diffpos_twobodyJ2) ' km']);
diffvel_twobodyJ2 = norm(validate_twobodyJ2(end, 4:6) - vf_esat);
disp(['Difference in Velovity = ' num2str(diffvel_twobodyJ2) ' km/s']);

% Integrate 3D two-body problem with J2, Sun-Earth, Earth-Moon perturbations
disp('Two Body Problem (J2, S-E, E-M)');
% Initial parameters (4/18/2022 16:00:00) (STK)
r0_esat = [40444.664764 -11913.011610 -59.307254]; % Earth-Satellite
v0_esat = [0.868736 2.949561 -0.002695];
X0 = [r0_esat'; v0_esat'; r0_se'; v0_se'; r0_em'; v0_em'];
% Final parameters (4/28/2022 16:00:00) (STK)
rf_esat = [41884.897207 -4831.335067 -75.522676];
vf_esat = [0.352306 3.054591 -0.001348];
% Propagate orbit
effects = [1 1 1];
[~, validate_twobodyJ2SM] = ode45(fn, tspan, X0, options);
% Print results
diffpos_twobodyJ2SM = norm(validate_twobodyJ2SM(end, 1:3) - rf_esat);
disp(['Difference in Position = ' num2str(diffpos_twobodyJ2SM) ' km']);
diffvel_twobodyJ2SM = norm(validate_twobodyJ2SM(end, 4:6) - vf_esat);
disp(['Difference in Velovity = ' num2str(diffvel_twobodyJ2SM) ' km/s']);

toc;
%% Validation - DirecTV-9S
tic;
disp('Validation - DirecTV-9S');

% Initial parameters (4/18/2022 16:00:00) (NASA JPL HORIZONS)
r0_se = [-1.322574708254664E+08 -7.120535511008461E+07 4.289022599551827E+03]; % Sun-Earth
v0_se = [1.362658131963260E+01 -2.633474893703986E+01 2.439047330110355E-03];
r0_em = [-2.177334626275463E+05 -2.938494594337748E+05 -6.810797786048497E+02]; % Earth-Moon
v0_em = [8.699532350602149E-01 -6.268144526653440E-01 -9.832388893488264E-02];

% Time span
tspan = 0:60:seconds(datetime(2022, 4, 28, 16, 0, 0) - datetime(2022, 4, 18, 16, 0, 0)); % sec

% Integrate 3D two-body problem
disp('Two Body Problem');

% Initial parameters (4/18/2022 16:00:00) (STK)
r0_esat = [40818.860981 -10534.330666 8.625958]; % Earth-Satellite
v0_esat = [0.767781 2.977876 -0.000200];
X0 = [r0_esat'; v0_esat'; r0_se'; v0_se'; r0_em'; v0_em'];
% Final parameters (4/28/2022 16:00:00) (STK)
rf_esat = [42012.968509 -3454.414083 8.035996];
vf_esat = [0.251422 3.065079 -0.000304];
% Propagate orbit
effects = [0 0 0];
[~, validate_twobody] = ode45(fn, tspan, X0, options);
% Print results
diffpos_twobody = norm(validate_twobody(end, 1:3) - rf_esat);
disp(['Difference in Position = ' num2str(diffpos_twobody) ' km']);
diffvel_twobody = norm(validate_twobody(end, 4:6) - vf_esat);
disp(['Difference in Velovity = ' num2str(diffvel_twobody) ' km/s']);

% Integrate 3D two-body problem with J2 perturbations
disp('Two Body Problem (J2)');
% Initial parameters (4/18/2022 16:00:00) (STK)
r0_esat = [40821.466003 -10524.226621 8.624941]; % Earth-Satellite
v0_esat = [0.767044 2.978066 -0.000200];
X0 = [r0_esat'; v0_esat'; r0_se'; v0_se'; r0_em'; v0_em'];
% Final parameters (4/28/2022 16:00:00) (STK)
rf_esat = [42029.470907 -3247.211536 8.004996];
vf_esat = [0.236307 3.066282 -0.000308];
% Propagate orbit
effects = [1 0 0];
[~, validate_twobodyJ2] = ode45(fn, tspan, X0, options);
% Print results
diffpos_twobodyJ2 = norm(validate_twobodyJ2(end, 1:3) - rf_esat);
disp(['Difference in Position = ' num2str(diffpos_twobodyJ2) ' km']);
diffvel_twobodyJ2 = norm(validate_twobodyJ2(end, 4:6) - vf_esat);
disp(['Difference in Velovity = ' num2str(diffvel_twobodyJ2) ' km/s']);

% Integrate 3D two-body problem with J2, Sun-Earth, Earth-Moon perturbations
disp('Two Body Problem (J2, S-E, E-M)');
% Initial parameters (4/18/2022 16:00:00) (STK)
r0_esat = [40764.167669 -10741.768485 -61.561203]; % Earth-Satellite
v0_esat = [0.782938 2.974008 -0.002342];
X0 = [r0_esat'; v0_esat'; r0_se'; v0_se'; r0_em'; v0_em'];
% Final parameters (4/28/2022 16:00:00) (STK)
rf_esat = [41999.078126 -3613.975495 -78.911345];
vf_esat = [0.263065 3.064176 -0.000996];
% Propagate orbit
effects = [1 1 1];
[~, validate_twobodyJ2SM] = ode45(fn, tspan, X0, options);
% Print results
diffpos_twobodyJ2SM = norm(validate_twobodyJ2SM(end, 1:3) - rf_esat);
disp(['Difference in Position = ' num2str(diffpos_twobodyJ2SM) ' km']);
diffvel_twobodyJ2SM = norm(validate_twobodyJ2SM(end, 4:6) - vf_esat);
disp(['Difference in Velovity = ' num2str(diffvel_twobodyJ2SM) ' km/s']);

toc;