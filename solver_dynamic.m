function solver_dynamic
% SOLVER_DYNAMIC solves the DAE and ODE of the Model for a variety of parameter combinations.
% The chosen dynamic contact angle model is taken from the molecular
% kinetics theory.
%
% The wall functions considered are:
%   * a 'constricted' case: w(x) = 2/3 + cos(2 pi x)/3
%   * a 'constant' case:    w(x) = 1
%
% (c) 2020 Stephan B. Lunowa
%
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

%% Parameters for constricted case
% final time (too big for most simulations, so they will stop early)
T = 2;
% maximal time step size
dt = 1e-2;

% wall function
w = @(x) 2.0/3.0 + cos(2 * pi * x) / 3.0;
% total flow
q = @(t) 0.6667 * ones(size(t));
% inlet pressure
p_in = @(t) 12 * ones(size(t));
% initial position of the interface
gamma0 = 0;

% folder to save solutions
folder = 'constricted/dynamic/';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process
eta = 0:0.25:1; % dynamic contact angle coeff.
%for i = 1:numel(eta)
parfor i = 1:numel(eta)
    for slip = [0, 1/6] % slip length
        for M = [0.5, 1, 2] % viscosity ratio (fluid I/II)
            for Ca = [0.5, 1, 2] % capillary number
                % contact angle function
                theta = @(x, u) acos(max(min(1/2 + eta(i) * Ca * u, 1), -1));

                % create model, solve and save
                m = Model(w, theta, Ca, M, slip);
                solve_and_save(m, p_in, q, gamma0, T, dt, folder, eta(i))
            end
        end
    end
end

%% Parameters for constant case
% final time (too big for most simulations, so they will stop early)
T = 2;
% maximal time step size
dt = 1e-2;

% wall function
w = @(x) ones(size(x));
% total flow
q = @(t) ones(size(t));
% inlet pressure
p_in = @(t) 3 * ones(size(t));
% initial position of the interface
gamma0 = 0;

% folder to save solutions
folder = 'constant/dynamic/';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process
eta = 0:0.25:1; % dynamic contact angle coeff.
%for i = 1:numel(eta)
parfor i = 1:numel(eta)
    for slip = [0, 1/6] % slip length
        for M = [0.5, 1, 2] % viscosity ratio (fluid I/II)
            for Ca = [0.5, 1, 2] % capillary number
                % contact angle function
                theta = @(x, u) acos(max(min(1/2 + eta(i) * Ca * u, 1), -1));

                % create model, solve and save
                m = Model(w, theta, Ca, M, slip);
                solve_and_save(m, p_in, q, gamma0, T, dt, folder, eta(i))
            end
        end
    end
end

end

function solve_and_save(m, p_in, q, gamma0, T, dtMax, folder, eta)
    % SOLVE_AND_SAVE helper function to automatize the solution and saving process

    % Solve dae model
    options = odeset('MaxStep', dtMax, 'RelTol',1e-3, 'AbsTol',[1e-6, 1e-7]); % solver options
    gammaMin = 1e-3;
    if m.M < 1e-2 && gamma0 < gammaMin
        m1 = m.solveDAE(p_in, T, gammaMin, options);
    else
        m1 = m.solveDAE(p_in, T, gamma0, options);
    end
    m1.saveSolution([folder '/dae_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '_eta' num2str(eta) '.dat'])

    % Solve ode model
    options = odeset('MaxStep', dtMax, 'RelTol',1e-3, 'AbsTol', 1e-6); % solver options
    m2 = m.solveODE(q, T, gamma0, options);
    m2.saveSolution([folder '/ode_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '_eta' num2str(eta) '.dat'])
end
