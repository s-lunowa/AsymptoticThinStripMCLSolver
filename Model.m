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
        % The computed saturation
        sat     {mustBeNumeric}
        % The computed local capillary pressure
        p_cloc  {mustBeNumeric}
        % The computed intrinsic pressure difference
        p_diff  {mustBeNumeric}
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
            int_fun = @(x) 1.0 ./ (m.w(x) .* m.w(x) .* (m.w(x) + 3 * m.slip));
            integral01 = integral(int_fun, 0, 1);
            dae = @(t, y, dydt) ...
                [p_in(t) - 3 * y(2) .* ((1 - m.M) * integral(int_fun, 0, y(1)) + m.M * integral01) ...
                    - cos(m.theta(y(1), y(2) ./ m.w(y(1)))) ./ (m.Ca * m.w(y(1))); ...
                 dydt(1) - y(2) ./ m.w(y(1))];
            % initial total flux (approx)
            q0 = (p_in(0) - cos(m.theta(gamma0, 0 / m.w(gamma0))) / (m.Ca * m.w(gamma0))) ...
               / (3 * m.M * integral01 + 3 * (1-m.M) * integral(int_fun, 0, gamma0));

            % solve model
            [y0, yp0] = decic(dae, 0, [gamma0, q0], [1 0], [q0/m.w(0) 0], [], options); % good initial conditions
            [m.time, sol] = ode15i(dae, [0, T], y0, yp0, options); % solve dae
            m.gamma = sol(:,1);
            m.q = sol(:,2);

            % compute and set additional data
            m.p_in = p_in(m.time);
            [m.sat, m.p_cloc, tau] = computeSatP_clocP_diff(m, m.gamma, m.q);
            m.p_diff = m.p_cloc + tau .* m.q / integral(m.w, 0,1);
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

            % solve model
            ode = @(t, y) q(t) ./ m.w(y);
            if nargin < 5 % solve ode w/ or w/o/ options given
                [m.time,m.gamma] = ode45(ode, [0,T], gamma0);
            else
                [m.time,m.gamma] = ode45(ode, [0,T], gamma0, options);
            end

            % compute and set additional data
            m.q = q(m.time);
            [m.sat, m.p_cloc, tau] = computeSatP_clocP_diff(m, m.gamma, m.q);
            m.p_diff = m.p_cloc + tau .* m.q / integral(m.w, 0,1);
            int_fun = @(x) 1.0 ./ (m.w(x) .* m.w(x) .* (m.w(x) + 3 * m.slip));
            integral01 = integral(int_fun, 0, 1);
            m.p_in = 3 * m.q .* ((1 - m.M) * arrayfun(@(x) integral(int_fun, 0,x), m.gamma) + m.M * integral01) ...
                   + m.p_cloc;
        end

        function [] = plot(m,dae)
            %PLOT plots the solution over time in one figure.
            %
            % INPUT:
            %   dae - whether the dae model was solved (bool, ode otherwise)

            fig_title = ['solution for M = ' num2str(m.M) ', slip = ' num2str(m.slip) ', Ca = ' num2str(m.Ca)];
            if dae
                figure('Name', ['DAE ' fig_title])
            else
                figure('Name', ['ODE ' fig_title])
            end

            stackedplot(m.time, [m.gamma, m.p_in, m.q, m.sat, m.p_cloc, m.p_diff], ...
                        'DisplayLabels',["gamma", "p_in", "q", "S", "p_cloc", "p_diff"])
        end

        function [] = plot_p_diff_sat(m)
            %PLOT plots the static p_c - sat and the (solution dependent) dynamic p_diff - sat curve in one figure.
            %

            figure('Name', ['p_diff - sat relation for M = ' num2str(m.M) ', slip = ' num2str(m.slip) ', Ca = ' num2str(m.Ca)])
            % dynamic p_diff - s
            plot(m.sat,m.p_diff)
            hold on
            % static p_c - s
            g = linspace(0,1, 100);
            % compute the saturation S = W(gamma) / W(1)
            W1 = integral(m.w, 0,1);
            s = arrayfun(@(g) integral(m.w, 0,g) ./ W1, g);
            % compute the local capillary pressure
            p_c = cos(m.theta(g, zeros(size(g)))) ./ (m.Ca * m.w(g));
            plot(s, p_c)
            hold off
            legend('dyn. p_{diff}', 'static p_c')
            xlabel('Saturation S')
            xlim([0,1])
        end

        function [] = saveSolution(m, filename)
            %SAVESOLUTION writes a csv table with the solution of the model.
            %
            % INPUT:
            %   filename - filename (default 'YYYY-MM-DD_HH-MM-SS.dat')

            % create filename if not given
            if nargin < 2
                t = datetime();
                t.Format = 'uuuu-MM-dd_HH-mm-ss';
                filename = string(t) + '.dat';
            end

            % write csv file
            fileID = fopen(filename,'w');
            fprintf(fileID,'# Ca = %g  M = %g  slip = %g  w = %s  theta = %s\n', ...
                    m.Ca, m.M, m.slip, func2str(m.w), strrep(func2str(m.theta), ',', ';'));
            fprintf(fileID, 'time,gamma,p_in,q,s,p_cloc,p_diff\n');
            fprintf(fileID, '%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
                            [m.time, m.gamma, m.p_in, m.q, m.sat, m.p_cloc, m.p_diff]');
            fclose(fileID);
        end

        function [] = saveCapillarityParameters(m, filename)
            %SAVECAPILLARITYPARAMETERS writes a csv table with the local static capillary pressure and the dynamic coefficient.
            %
            % INPUT:
            %   filename - filename (default 'YYYY-MM-DD_HH-MM-SS.dat')

            % compute values
            g = linspace(0,1, 200)';
            [s, p_c, tau] = computeSatP_clocP_diff(m, g, zeros(size(g)));

            % create filename if not given
            if nargin < 2
                t = datetime();
                t.Format = 'uuuu-MM-dd_HH-mm-ss';
                filename = string(t) + '.dat';
            end

            % write csv file
            fileID = fopen(filename,'w');
            fprintf(fileID,'# Ca = %g  M = %g  slip = %g  w = %s  theta = %s\n', ...
                    m.Ca, m.M, m.slip, func2str(m.w), strrep(func2str(m.theta), ',', ';'));
            fprintf(fileID,'s,p_cloc,tau\n');
            fprintf(fileID,'%.6g,%.6g,%.6g\n', [s, p_c, tau]');
            fclose(fileID);
        end
    end

    methods (Access = protected)
        function [sat, p_cloc, tau] = computeSatP_clocP_diff(m, gamma, q)
            %COMPUTESATP_CLOCP_TAU computes the saturation, the local capillary pressure and the dynamic coefficient from a solution gamma, p_in, q.
            %
            % The input may be vectorial, but must have the same size.
            %
            % INPUT:
            %   gamma - the interface position
            %   q     - the total flux

            % compute the saturation S = W(gamma) / W(1)
            W1 = integral(m.w, 0,1);
            sat = arrayfun(@(g) integral(m.w, 0,g) ./ W1, gamma);

            % compute the local capillary pressure
            wg = m.w(gamma);
            p_cloc = cos(m.theta(gamma, q ./ wg)) ./ (m.Ca * wg);

            % compute the dynamic coefficient
            % tau = 3 W(1) (\int_0^gamma W(x) / (w(x))^2 (3 slip + w(x))) dx / W(gamma)
            %              + m.M * \int_gamma^1 (W(1) - W(x)) / ((w(x))^2 (3 slip + w(x))) dx / (W(1) - W(gamma)))
            tmp1 = @(x) integral(m.w, 0,x) ./ ((m.w(x)).^2 .* (3 * m.slip + m.w(x)));
            tmp2 = @(x) integral(m.w, x,1) ./ ((m.w(x)).^2 .* (3 * m.slip + m.w(x)));
            int1 = @(g) integral(@(x) arrayfun(tmp1,x), 0,g) ./ integral(m.w, 0,g);
            int2 = @(g) integral(@(x) arrayfun(tmp2,x), g,1) ./ integral(m.w, g,1);
            tau = arrayfun(@(g) 3 * W1 * (int1(g) + m.M * int2(g)), gamma);
            % correct for gamma == 0 or gamma == 1 (int1 == 0 or int2 == 0)
            tau(gamma == 0) = 3 * W1 * m.M * int2(0);
            tau(gamma == 1) = 3 * W1 * int1(1);
        end
    end
end
