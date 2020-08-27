function solver_constant
% SOLVER_CONSTANT solves the DAE and ODE of the Model for a variety of parameter combinations.
% The chosen contact angle model is static.
%
% The wall function considered is constant w(x) = 1
%
% (c) 2020 Stephan B. Lunowa
%
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.


%% Parameters
% final time (too big for most simulations, so they will stop early)
T = 1.5;
% maximal time step size
dtMax = 1e-2;

% total flow
q = @(t) ones(size(t));
% wall function
w = @(x) ones(size(x));
% contact angle function
theta = @(x, u) ones(size(x)) * (pi / 3);
% initial position of the interface
gamma0 = 0;

% folder to save solutions
folder = 'constant';
if ~exist(folder, 'dir')
   mkdir(folder)
end

%% Solution process

M = [0, 1e-3, 0.1, 0.5, 1, 2]; % viscosity ratio (fluid I/II)
%for i = 1:numel(M)
parfor i = 1:numel(M)
    for slip = [0, 0.05, 1/6, 0.5] % slip length
        for Ca = [0.1, 0.25, 0.5, 1, 2] % capillary number
            % inlet pressure
            p_tmp = 2 + cos(theta(0,0)) / Ca;
            p_in = @(t) p_tmp * ones(size(t));

            % create, solve and save model
            m = Model(w, theta, Ca, M(i), slip);
            solve_and_save(m, p_in, q, gamma0, T, dtMax, folder)
        end
    end
end
end

function solve_and_save(m, p_in, q, gamma0, T, dtMax, folder)
    % SOLVE_AND_SAVE helper function to automatize the solution and saving process

    m.saveCapillarityParameters([folder '/s-p_c-tau_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

    % Solve dae model
    options = odeset('MaxStep', dtMax, 'RelTol',1e-3, 'AbsTol',[1e-6, 1e-7]); % solver options
    gammaMin = 1e-3;
    if m.M < 1e-2 && gamma0 < gammaMin
        m1 = m.solveDAE(p_in, T, gammaMin, options);
    else
        m1 = m.solveDAE(p_in, T, gamma0, options);
    end
    m1.saveSolution([folder '/dae_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])

    % Solve ode model
    options = odeset('MaxStep', dtMax, 'RelTol',1e-3, 'AbsTol', 1e-6); % solver options
    m2 = m.solveODE(q, T, gamma0, options);
    m2.saveSolution([folder '/ode_M' num2str(m.M) '_slip' num2str(m.slip) '_Ca' num2str(m.Ca) '.dat'])
end
