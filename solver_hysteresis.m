function solver_hysteresis
% SOLVER_HYSTERESIS solves the DAE and ODE of the Model for a variety of parameter combinations.
% The chosen contact angle models are hysteretic, dynamic and static.
%
% The wall functions considered are:
%   * a 'constant' case:    w(x) = 1
%   * a 'constricted' case: w(x) = 2/3 + cos(2 pi x)/3
%
% (c) 2021 Stephan B. Lunowa
%
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

%% Parameters for the constant width case
% final time (too big for most simulations, so they will stop early)
T = 4.5;
% maximal time step size
dt = 1e-2;

% wall function
w = @(x) ones(size(x));
% slip length
slip = 1/6;
% viscosity ratio (fluid I/II)
M = 1;
% capillary number
Ca = 0.5;
% static contact angle
theta_static = pi / 3;
% deviation of hysteretic to static contact angle
theta_hyst = pi / 12;

% initial position of the interface
gamma0 = 0;

% inlet pressure
p_in = @(t) cos(theta_static)/Ca + 2 * (1 - t/2); % constant

% folder to save solutions
folder = 'constant/hysteresis/';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process
eta = 0:0.5:1; % dynamic contact angle coeff.
for i = 1:numel(eta)
%parfor i = 1:numel(eta)
    % contact angle function
    c1 = cos(theta_static);
    c2 = cos(theta_static + theta_hyst);
    c3 = cos(theta_static - theta_hyst);
    theta = @(x, u) acos(max(min(c1 + eta(i) * Ca * u, 1), -1));
    zeta = @(x, p) min(p - c2, max(0, p - c3)) / (eta(i) * Ca);

    % create model, solve and save
    m = Model(w, theta, Ca, M, slip);
    m_hyst = ModelHysteretic(w, zeta, Ca, M, slip);
    solve_and_save(m, m_hyst, p_in, gamma0, T, dt, folder, eta(i))
end

%% Parameters for the constricted width case
% final time (too big for most simulations, so they will stop early)
T = 4.5;
% maximal time step size
dt = 1e-2;

% wall function
w = @(x) 2/3.0 + cos(2 * pi * x) / 3.0;
% slip length
slip = 1/6;
% viscosity ratio (fluid I/II)
M = 1;
% capillary number
Ca = 0.5;
% static contact angle
theta_static = pi / 3;
% deviation of hysteretic to static contact angle
theta_hyst = pi / 12;

% initial position of the interface
gamma0 = 0;

% inlet pressure
%p_in = @(t) cos(theta_static)/Ca + 6.25 * cos(pi/4 * t); % constant
p_in = @(t) cos(theta_static)/Ca + 8 * (1 - t/2); % constant

% folder to save solutions
folder = 'constricted/hysteresis/';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process
eta = 0:0.5:1; % dynamic contact angle coeff.
for i = 1:numel(eta)
%parfor i = 1:numel(eta)
    % contact angle function
    c1 = cos(theta_static);
    c2 = cos(theta_static + theta_hyst);
    c3 = cos(theta_static - theta_hyst);
    theta = @(x, u) acos(max(min(c1 + eta(i) * Ca * u, 1), -1));
    zeta = @(x, p) min(p - c2, max(0, p - c3)) / (eta(i) * Ca);

    % create model, solve and save
    m = Model(w, theta, Ca, M, slip);
    m_hyst = ModelHysteretic(w, zeta, Ca, M, slip);
    solve_and_save(m, m_hyst, p_in, gamma0, T, dt, folder, eta(i))
end

end

function solve_and_save(m, m_hyst, p_in, gamma0, T, dtMax, folder, eta)
    % SOLVE_AND_SAVE helper function to automatize the solution and saving process

    % Solve dae model
    options = odeset('MaxStep', dtMax, 'RelTol',1e-3, 'AbsTol',[1e-6, 1e-7]); % solver options
    gammaMin = 1e-3;
    if m.M < 1e-2 && gamma0 < gammaMin
        gamma0 = gammaMin;
    end
    m1 = m.solveDAE(p_in, T, gamma0, options);
    m1.saveSolution([folder '/dae_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '_eta' num2str(eta) '.dat'])

    if eta > 0
        m2 = m_hyst.solveDAE_explicit(p_in, T, gamma0, ceil(T/dtMax));
        m2.saveSolution([folder '/dae_hyst_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '_eta' num2str(eta) '.dat'])
    end
end