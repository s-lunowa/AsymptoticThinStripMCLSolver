classdef Model
%Model Contains all the data necessary for the solution of the thin-strip model.
%   This model can solve the dae/ode model
%
%       p_in(t) - 3 q(t) ((1 - M) int_0^gamma(t) ... dx + M * int_0^1 ... dx)
%           == cos(theta(gamma(t), q(t)/w(gamma(t)))) / (Ca w(gamma(t)))
%       dgamma(t)/dt == q(t) / w(gamma(t))
%
% over time [0,T] with respect to gamma(0) == gamma0.
%
% (c) 2020 Stephan B. Lunowa
%
% This software is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

    properties
        % Capillary number
        Ca  (1,1) {mustBeNumeric,mustBePositive   ,mustBeFinite} = 1
        % Viscosity ratio (fluid I/II)
        M   (1,1) {mustBeNumeric,mustBeNonnegative,mustBeFinite} = 1
        % Slip length
        slip(1,1) {mustBeNumeric,mustBeNonnegative,mustBeFinite} = 0
        % Wall function (depends on x, uniformly positive) [vectorized]
        w    (1,1) function_handle = @(x) ones(size(x))
        % Contact angle function (depends on x and u) [vectorized]
        theta(1,1) function_handle = @(x,u) pi/2 * ones(size(x))

        % The solution for the interface position
        gamma   {mustBeNumeric}
        % The solution for the inlet pressure
        p_in    {mustBeNumeric}
        % The solution for the total flux
        q       {mustBeNumeric}
        % The simulated time points
        time    {mustBeNumeric}
    end

    methods
        function m = Model(w, theta, Ca, M, slip)
            %MODEL Constructs an instance of the model.
            %
            % INPUT:
            %   w     - wall function (depends on x, uniformly positive, default constant 1) [vectorized]
            %   theta - contact angle function (depends on x and u, default constant pi/2) [vectorized]
            %   Ca    - capillary number (positive, default 1)
            %   M     - viscosity ratio (fluid I/II, non-negative, default 1)
            %   slip  - slip coefficient (non-negative, default 0)

            if nargin > 0
                if abs(w(0)-1) > 1e-6
                    error('The wall function w is not 1 at x=0.')
                end
                m.w = w;
                if nargin > 1
                    m.theta = theta;
                    if nargin > 2
                        m.Ca = Ca;
                        if nargin > 3
                            m.M = M;
                            if nargin > 4
                                m.slip = slip;
                            end
                        end
                    end
                end
            end
        end

        function m = solveDAE(m, p_in, T, gamma0, options)
            % SOLVEDAE solves the dae model. The solution is in the properties of the returned Model.
            %
            % INPUT:
            %   p_in    - inlet pressure function (depends on t)
            %   T       - final time of simulation time [0,T] (positive)
            %   gamma0  - initial interface position (non-negative)
            %   options - odeset for dae-solver options (otional)

            if T <= 0.0
                error('The final time T must be positive.')
            end
            if (gamma0 < 0.0 || gamma0 > 1.0)
                error('The initial interface position gamma0 must be non-negative.')
            end
            if nargin < 5
                options = odeset();
            else
                options = odeset(options, 'Jacobian', {[],[0 0; 1 0]});
            end

            % create model
            integral_fun = @(x) 1.0 ./ (m.w(x) .* m.w(x) .* (m.w(x) + 3 * m.slip));
            integral01 = integral(integral_fun, 0, 1);
            dae = @(t, y, dydt) ...
                [p_in(t) - 3 * y(2) .* ((1 - m.M) * integral(integral_fun, 0, y(1)) + m.M * integral01) ...
                    - cos(m.theta(y(1), y(2) ./ m.w(y(1)))) ./ (m.Ca * m.w(y(1))); ...
                 dydt(1) - y(2) ./ m.w(y(1))];
            % initial total flux (approx)
            q0 = (p_in(0) - cos(m.theta(gamma0, 0 / m.w(gamma0))) / (m.Ca * m.w(gamma0))) ...
               / (3 * m.M * integral01 + 3 * (1-m.M) * integral(integral_fun, 0, gamma0));

            % solve model
            [y0, yp0] = decic(dae, 0, [gamma0, q0], [1 0], [q0/m.w(0) 0], [], options); % good initial conditions
            [m.time, sol] = ode15i(dae, [0, T], y0, yp0, options); % solve dae
            m.gamma = sol(:,1);
            m.q = sol(:,2);
            m.p_in = p_in(m.time);
        end

        function m = solveODE(m, q, T, gamma0, options)
            % SOLVEODE solves the ode model. The solution is in the properties of the returned Model.
            %
            % INPUT:
            %   q       - total flux function (depends on t)
            %   T       - final time of simulation time [0,T] (positive)
            %   gamma0  - initial interface position (non-negative)
            %   options - odeset for dae-solver options (otional)

            if T <= 0.0
                error('The final time T must be positive.')
            end
            if (gamma0 < 0.0 || gamma0 > 1.0)
                error('The initial interface position gamma0 must be non-negative.')
            end

            % create model
            integral_fun = @(x) 1.0 ./ (m.w(x) .* m.w(x) .* (m.w(x) + 3 * m.slip));
            integral01 = integral(integral_fun, 0, 1);

            % solve model
            ode = @(t, y) q(t) ./ m.w(y);
            if nargin < 5 % solve ode w/ or w/o/ options given
                [m.time,m.gamma] = ode45(ode, [0,T], gamma0);
            else
                [m.time,m.gamma] = ode45(ode, [0,T], gamma0, options);
            end
            m.q = q(m.time);
            m.p_in = cos(m.theta(m.gamma, m.q ./ m.w(m.gamma))) ./ (m.Ca * m.w(m.gamma));
            for i = 1:length(m.time)
                m.p_in(i) = m.p_in(i) + 3 * m.q(i) * ((1 - m.M) * integral(integral_fun, 0,m.gamma(i)) + m.M * integral01);
            end
        end

        function [] = plot(m,dae)
            %PLOT plots the solution over time in one figure.
            %
            % INPUT:
            %   dae - whether the dae model was solved (bool, ode otherwise)

            if dae
                figure('Name', 'Solution for dae model')
            else
                figure('Name', 'Solution for ode model')
            end

            yyaxis left
            plot(m.time,m.gamma)
            xlabel('time t')
            ylabel('\gamma(t)')

            yyaxis right
            if(dae)
                plot(m.time, m.q)
                ylabel('q(t)')
            else
                plot(m.time, m.p_in)
                ylabel('p_{in}(t)')
            end
        end

        function [] = saveSolution(m, filename)
            %SAVESOLUTION writes a csv table with the solution of the model.
            %
            % INPUT:
            %   filename - filename (default 'YYYY-MM-DD_HH-MM-SS.dat')

            time = m.time; gamma = m.gamma; p_in = m.p_in; q = m.q;
            param = {'Ca'; num2str(m.Ca); 'M'; num2str(m.M); ...
                     'slip'; num2str(m.slip); 'w'; func2str(m.w); ...
                     'theta'; strrep(func2str(m.theta), ',', ';')};
            for i = 11:length(time)
                 param{i} = ' ';
            end
            T = table(time, gamma, p_in, q, param);
            if nargin < 2
                t = datetime();
                t.Format = 'uuuu-MM-dd_HH-mm-ss';
                filename = string(t) + '.dat';
            end
            writetable(T, filename)
        end

        function [] = saveCapillarityParameters(m, filename)
            %SAVECAPILLARITYPARAMETERS writes a csv table with the local static capillary pressure and the dynamic coefficient.
            %
            % INPUT:
            %   filename - filename (default 'YYYY-MM-DD_HH-MM-SS.dat')

            tmp1 = @(a) integral(m.w, 0,a) ./ ((m.w(a)).^2 .* (3 * m.slip + m.w(a)));
            tmp2 = @(a) integral(m.w, a,1) ./ ((m.w(a)).^2 .* (3 * m.slip + m.w(a)));
            integ1 = @(a) integral(tmp1, 0,a, 'ArrayValued', true) ./ integral(m.w, 0,a);
            integ2 = @(a) integral(tmp2, a,1, 'ArrayValued', true) ./ integral(m.w, a,1);
            W1 = integral(m.w, 0,1);

            tau_g = @(g) 3 * W1 * (integ1(g) + m.M * integ2(g));
            Psi_g = @(g) 1 ./ m.w(g);
            sat = @(g) integral(m.w, 0,g) ./ W1;

            g = linspace(0,1, 200)';
            s = zeros(size(g));
            tau = zeros(size(g));
            p_cloc = zeros(size(g));
            parfor i = 1:length(g)
                s(i) = sat(g(i));
                tau(i) = tau_g(g(i));
                p_cloc(i) = cos(m.theta(g(i), 0)) / m.Ca * Psi_g(g(i));
            end
            param = {'Ca'; num2str(m.Ca); 'M'; num2str(m.M); ...
                     'slip'; num2str(m.slip); 'w'; func2str(m.w); ...
                     'theta'; strrep(func2str(m.theta), ',', ';')};
            for i = 11:length(g)
                 param{i} = ' ';
            end
            T = table(s, p_cloc, tau, param);
            if nargin < 2
                t = datetime();
                t.Format = 'uuuu-MM-dd_HH-mm-ss';
                filename = string(t) + '.dat';
            end
            writetable(T, filename)
        end
    end
end
