%%
% (c) 2020 Stephan B. Lunowa
%
% This software is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

function solver_constricted

%% Parameters for constricted case
T = 0.6667; % final time
NT = 100;   % minimal number of time steps

% total flow
q = @(t) ones(size(t));
% wall function
w = @(x) 2.0/3.0 + cos(2 * pi * x) / 3.0;
% contact angle function
theta = @(x, u) ones(size(x)) * (pi / 3);
% initial position of the interface
gamma0 = 0;

% folder to save solutions
folder = 'constricted';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process
for slip = [0, 0.025, 0.05, 1/6, 0.5] % slip length
    for M = [0, 1e-6, 0.1, 0.5, 1, 2] % viscosity ratio (fluid I/II)
        for Ca = [0.1, 0.25, 0.5, 1, 2] % capillary number
            % inlet pressure
            p_tmp = 6 * (M + 1) / (slip + 1) + 1 / Ca;
            p_in = @(t) p_tmp * ones(size(t));

            % create model
            m = Model(w, theta, Ca, M, slip);
            m.saveCapillarityParameters([folder '/s-p_c-tau_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

            % Solve and save
            solve_and_save(m, p_in, q, gamma0, T, T/NT, folder)
        end
     	%pause
    end
    close all
end

%% Parameters for long constricted case
T = 0.5;      % final time
NT = 100;   % minimal number of time steps

% inlet pressure
p_in = @(t) 10 * ones(size(t)); % 3 + 2 * t;
% total flow
q = @(t) ones(size(t));
% wall function
w = @(x) 2.0/3.0 + (cos(4 * pi * x) .* (x <= 0.25 | x >= 0.75) - (0.25 < x) .* (x < 0.75)) / 3.0;
% contact angle function
theta = @(x, u) ones(size(x)) * (pi / 3);

% initial position of the interface
gamma0 = 0;

% folder to save solutions
folder = 'constrictedLong';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process
for slip = [0, 0.025, 0.05, 1/6, 0.5] % slip length
    for M = [0, 1e-6, 0.1, 0.5, 1, 2] % viscosity ratio (fluid I/II)
        for Ca = [0.1, 0.25, 0.5, 1, 2] % capillary number
            % inlet pressure
            p_tmp = 12 * (M + 1) / (slip + 1) + 1 / Ca;
            p_in = @(t) p_tmp * ones(size(t));

            % create model
            m = Model(w, theta, Ca, M, slip);
            m.saveCapillarityParameters([folder '/s-p_c-tau_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

            % Solve and save
            solve_and_save(m, p_in, q, gamma0, T, T/NT, folder)
        end
     	%pause
    end
    close all
end
end

% helper function to automatize the solution and saving process
function solve_and_save(m, p_in, q, gamma0, T, MaxStep, folder)
    % Solve dae model
    options = odeset('MaxStep', MaxStep, 'RelTol',1e-3, 'AbsTol',[1e-6, 1e-7]); % solver options
    if m.M < 1e-6 && gamma0 < 5e-2
        m1 = m.solveDAE(p_in, T, 5e-2, options);
    elseif gamma0 < 1e-3
        m1 = m.solveDAE(p_in, T, 1e-3, options);
    else
        m1 = m.solveDAE(p_in, T, gamma0, options);
    end
    m1.plot(true)
    title(['DAE with M = ' num2str(m.M) ', \lambda = ' num2str(m.slip) ', Ca = ' num2str(m.Ca)])
    m1.saveSolution([folder '/dae_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

    % Solve ode model
    options = odeset('MaxStep', MaxStep, 'RelTol',1e-3, 'AbsTol', 1e-6); % solver options
    m2 = m.solveODE(q, T, gamma0, options);
    m2.plot(false)
    title(['ODE with M = ' num2str(m.M) ', \lambda = ' num2str(m.slip) ', Ca = ' num2str(m.Ca)])
    m2.saveSolution([folder '/ode_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])
end
