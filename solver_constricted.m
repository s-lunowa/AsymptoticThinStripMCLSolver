function solver_constricted
% SOLVER_CONSTRICTED solves the DAE and ODE of the Model for a variety of parameter combinations.
% The chosen contact angle model is static.
%
% The wall functions considered are:
%   * a 'constricted' case:      w(x) = 2/3 + cos(2 pi x)/3
%   * a 'long constricted' case: w(x) = 2/3 + (cos(4 * pi * x)/3 if (x <= 0.25 or x >= 0.75)
%                                     = 1/3 otherwise;
%
% (c) 2020 Stephan B. Lunowa
%
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.


%% Parameters for constricted case
% final time (too big for most simulations, so they will stop early)
T = 1.5;
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
for M = [0, 1e-6, 0.1, 0.5, 1, 2] % viscosity ratio (fluid I/II)
    for slip = [0, 0.025, 0.05, 1/6, 0.5] % slip length
        disp(['M = ' num2str(M) ', slip = ' num2str(slip)])
        % capillary number
        Ca = [0.1, 0.25, 0.5, 1, 2];
        parfor i = 1:length(Ca)
        %for i = 1:length(Ca)
            % inlet pressure
            p_tmp = 6 * (M + 1) / (slip + 1) + cos(theta(0,0)) / Ca(i);
            p_in = @(t) p_tmp * ones(size(t));

            % create, solve and save model
            m = Model(w, theta, Ca(i), M, slip);
            solve_and_save(m, p_in, q, gamma0, T, T/NT, folder)
        end
    end
end

%% Parameters for long constricted case
% final time (too big for most simulations, so they will stop early)
T = 1;
NT = 100;   % minimal number of time steps

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
for M = [0, 1e-6, 0.1, 0.5, 1, 2] % viscosity ratio (fluid I/II)
    for slip = [0, 0.025, 0.05, 1/6, 0.5] % slip length
        disp(['M = ' num2str(M) ', slip = ' num2str(slip)])
        parfor Ca = [0.1, 0.25, 0.5, 1, 2] % capillary number
            % inlet pressure
            p_tmp = 12 * (M + 1) / (slip + 1) + cos(theta(0,0)) / Ca;
            p_in = @(t) p_tmp * ones(size(t));

            % create, solve and save model
            m = Model(w, theta, Ca, M, slip);
            solve_and_save(m, p_in, q, gamma0, T, T/NT, folder)
        end
    end
end
end

function solve_and_save(m, p_in, q, gamma0, T, MaxStep, folder)
    % SOLVE_AND_SAVE helper function to automatize the solution and saving process

    m.saveCapillarityParameters([folder '/s-p_c-tau_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

    % Solve dae model
    options = odeset('MaxStep', MaxStep, 'RelTol',1e-3, 'AbsTol',[1e-6, 1e-7]); % solver options
    gammaMin = max(min(-log10(m.M),10),0.1) * 1e-2;
    if m.M < 1 && gamma0 < gammaMin
        m1 = m.solveDAE(p_in, T, gammaMin, options);
    else
        m1 = m.solveDAE(p_in, T, gamma0, options);
    end
    m1.saveSolution([folder '/dae_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

    % Solve ode model
    options = odeset('MaxStep', MaxStep, 'RelTol',1e-3, 'AbsTol', 1e-6); % solver options
    m2 = m.solveODE(q, T, gamma0, options);
    m2.saveSolution([folder '/ode_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])
end
