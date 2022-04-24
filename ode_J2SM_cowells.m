function Xdot = ode_J2SM_cowells(~, X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function integrates orbits using Cowell propagation
% and J2, Sun-Earth, and Earth-Moon perturbations
%
% Inputs
%   X - (18 x 1) Array of initial state vectors (x, y, z, vx, vy, vz) for a
%       satellite, the Sun, and the Moon
%
% Outputs
%   Xdot - (18 x 1) Array of time derivatives of state vectors 
%          (x, y, z, vx, vy, vz) considering J2, Sun-Earth, and Earth-Moon
%          perturbations
%
% Example
%
%   This functions is used within the built-in ode45 function 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare global variables
global mu Re J2 mu_s mu_m effects;

x_esat = X(1); % Earth-Satellite
y_esat = X(2);
z_esat = X(3);
vx_esat = X(4);
vy_esat = X(5);
vz_esat = X(6);
x_se = X(7); % Sun-Earth
y_se = X(8);
z_se = X(9);
vx_se = X(10);
vy_se = X(11);
vz_se = X(12);
x_em = X(13); % Earth-Moon
y_em = X(14);
z_em = X(15);
vx_em = X(16);
vy_em = X(17);
vz_em = X(18);

% Position vectors (km)
r_esat = sqrt(x_esat^2 + y_esat^2 + z_esat^2); % Earth-Satellite
r_se = sqrt(x_se^2 + y_se^2 + z_se^2); % Sun-Earth
r_me = sqrt((-x_em)^2 + (-y_em)^2 + (-z_em)^2); % Moon-Earth

if length(effects) ~= 3
    error(['Not accounting for the right amount of perturbation effects. ' ...
        'This function only calculates J2, Sun-Earth, and Earth-Moon effects.']);
end

xEffects = 0;
yEffects = 0;
zEffects = 0;

% J2 effects
if effects(1)
    J2x = (3*mu/(2*r_esat^5))*Re^2*J2*(1 - 5*(z_esat/r_esat)^2)*x_esat;
    J2y = (3*mu/(2*r_esat^5))*Re^2*J2*(1 - 5*(z_esat/r_esat)^2)*y_esat;
    J2z = (3*mu/(2*r_esat^5))*Re^2*J2*(3 - 5*(z_esat/r_esat)^2)*z_esat;
    
    xEffects = xEffects + J2x;
    yEffects = yEffects + J2y;
    zEffects = zEffects + J2z;
end

% Sun perturbation effects
if effects(2)
    r_ssat = [x_se y_se z_se] + [x_esat y_esat z_esat];
    x_ssat = r_ssat(1);
    y_ssat = r_ssat(2);
    z_ssat = r_ssat(3);
    SunEffectx = mu_s*((x_ssat/norm(r_ssat)^3) - (x_se/r_se^3));
    SunEffecty = mu_s*((y_ssat/norm(r_ssat)^3) - (y_se/r_se^3));
    SunEffectz = mu_s*((z_ssat/norm(r_ssat)^3) - (z_se/r_se^3));
    
    xEffects = xEffects + SunEffectx;
    yEffects = yEffects + SunEffecty;
    zEffects = zEffects + SunEffectz;
end

% Moon perturbation effects
if effects(3)
    r_msat = -[x_em y_em z_em] + [x_esat y_esat z_esat];
    x_msat = r_msat(1);
    y_msat = r_msat(2);
    z_msat = r_msat(3);
    MoonEffectx = mu_m*((x_msat/norm(r_msat)^3) - (-x_em/r_me^3));
    MoonEffecty = mu_m*((y_msat/norm(r_msat)^3) - (-y_em/r_me^3));
    MoonEffectz = mu_m*((z_msat/norm(r_msat)^3) - (-z_em/r_me^3));
    
    xEffects = xEffects + MoonEffectx;
    yEffects = yEffects + MoonEffecty;
    zEffects = zEffects + MoonEffectz;
end

% TODO: Solar radiation pressure
% if effects(4)
%     solradpress = 4.5e-6; % N/m^2
%     Ixx = 4500; % kg*m^2 (STK)
%     m = 2161; % kg (Vespucci paper)
%     r = sqrt(Ixx/(m*(2/5))); % m (model sat as solid sphere)
%     gamma = 1; % reflectivity (1 = sat absorbs all incident radiation)
%     
% end

Xdot = zeros(18, 1);

Xdot(1, 1) = vx_esat; % Earth-Satellite
Xdot(2, 1) = vy_esat;
Xdot(3, 1) = vz_esat;
Xdot(4, 1) = -mu*x_esat/r_esat^3 - xEffects;
Xdot(5, 1) = -mu*y_esat/r_esat^3 - yEffects;
Xdot(6, 1) = -mu*z_esat/r_esat^3 - zEffects;
Xdot(7, 1) = vx_se; % Sun-Earth
Xdot(8, 1) = vy_se;
Xdot(9, 1) = vz_se;
Xdot(10, 1) = -mu_s*x_se/r_se^3;
Xdot(11, 1) = -mu_s*y_se/r_se^3;
Xdot(12, 1) = -mu_s*z_se/r_se^3;
Xdot(13, 1) = vx_em; % Earth-Moon
Xdot(14, 1) = vy_em;
Xdot(15, 1) = vz_em;
Xdot(16, 1) = -mu_m*x_em/r_me^3; % r_me = r_em
Xdot(17, 1) = -mu_m*y_em/r_me^3;
Xdot(18, 1) = -mu_m*z_em/r_me^3;
end



%% BACKUP

% global mu1 mu2 m1 m2 R n;
% 
% x = X(1);
% y = X(2);
% z = X(3);
% vx = X(4);
% vy = X(5);
% vz = X(6);
% 
% Xdot = zeros(6, 1);
% 
% Xdot(1, 1) = vx;
% Xdot(2, 1) = vy;
% Xdot(3, 1) = vz;
% Xdot(4, 1) = n^2*x ...
%     + (-mu1*(x + (m2/(m1 + m2))*R)*((x + (m2/(m1 + m2))*R)^2 + y^2 + z^2)^(-3/2) ...
%     - mu2*(x - (m1/(m1 + m2))*R)*((x - (m1/(m1 + m2))*R)^2 + y^2 + z^2)^(-3/2)) ...
%     + 2*n*vy;
% Xdot(5, 1) = n^2*y ...
%     + (-mu1*(y)*((x + (m2/(m1 + m2))*R)^2 + y^2 + z^2)^(-3/2) ...
%     - mu2*(y)*((x - (m1/(m1 + m2))*R)^2 + y^2 + z^2)^(-3/2)) ...
%     - 2*n*vx;
% Xdot(6, 1) = -mu1*(z)*((x + (m2/(m1 + m2))*R)^2 + y^2 + z^2)^(-3/2) ...
%     - mu2*(z)*((x - (m1/(m1 + m2))*R)^2 + y^2 + z^2)^(-3/2);
% end