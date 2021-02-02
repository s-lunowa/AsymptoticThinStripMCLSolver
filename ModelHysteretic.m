classdef ModelHysteretic
%MODELHYSTERETIC contains all the data necessary for the solution of
% the thin-strip model involving a hysteretic contact angle model.
%
%   This model can solve the dae/ode model
%
%     zeta(gamma(t), w(gamma(t)) Ca (p_in(t) - q(t) J(gamma(t)))) == q(t) / w(gamma(t))
%     dgamma(t)/dt == q(t) / w(gamma(t))
%
%   where
%
%     J(g) = int_0^1 3/(w(x)^2 (w(x) + 3 slip)) dx + (M-1) int_g^1 3/(w(x)^2 (w(x) + 3 slip)) dx
%
%   over time [0,T] with respect to gamma(0) == gamma0.
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
        zeta(1,1) function_handle = @(x,p) min(p + 1, max(0, p - 1))

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
        function m = ModelHysteretic(w, zeta, Ca, M, slip)
            % MODELHYSTERETIC Constructs an instance of the model.
            %
            % INPUT:
            %   w     - wall function (depends on x, uniformly positive, default constant 1) [vectorized]
            %   zeta  - the inverse contact angle function (depends on x and p, default contact angle almost constant pi/2) [vectorized]
            %   Ca    - capillary number (positive, default 1)
            %   M     - viscosity ratio (fluid I/II, non-negative, default 1)
            %   slip  - slip coefficient (non-negative, default 0)

            if nargin > 0
                if abs(w(0)-1) > 1e-6
                    error('The wall function w is not 1 at x=0.')
                end
                m.w = w;
                if nargin > 1
                    m.zeta = zeta;
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
            % SOLVEDAE solves the dae model. The solution is in the properties of the returned ModelHysteretic.
            %
            % The solution process uses the system for dgamma/dt and q.
            %
            % INPUT:
            %   p_in    - inlet pressure function (depends on t)
            %   T       - final time of simulation time [0,T] (positive)
            %   gamma0  - initial interface position (0 <= gamma0 <= 1)
            %   options - odeset for dae-solver options (otional)

            if T <= 0.0
                error('The final time T must be positive.')
            end
            if (gamma0 < 0.0 || gamma0 > 1.0)
                error('The initial interface position gamma0 must be between 0 and 1.')
            end
            if nargin < 5
                options = odeset();
            end
            options = odeset(options, 'Jacobian', {[],[1 0; 0 0]}, 'Events', @ModelHysteretic.stopAtEnds);

            % create model
            int_fun = @(x) 3.0 ./ (m.w(x).^2 .* (m.w(x) + 3 * m.slip));
            integral01 = integral(int_fun, 0, 1);
            dae = @(t, y, dydt) ...
                [dydt(1) - y(2) ./ m.w(y(1)); ...
                 m.zeta(y(1), m.w(y(1)) * m.Ca * (p_in(t) - y(2) .* (integral01 + (m.M - 1) * integral(int_fun, y(1), 1)))) ...
                    - y(2) ./ m.w(y(1))
                ];
            % initial total flux (approx)
            tmp = integral01 + (m.M - 1) * integral(int_fun, gamma0,1);
            q0 = fzero(@(q) m.zeta(gamma0, m.w(gamma0) * m.Ca * (p_in(0) - q .* tmp)) - q ./ m.w(gamma0), 0);

            % solve model
            [y0, yp0] = decic(dae, 0, [gamma0, q0], [1 0], [q0/m.w(gamma0) 0], [], options); % good initial conditions
            [m.time, sol] = ode15i(dae, [0, T], y0, yp0, options); % solve dae
            m.gamma = sol(:,1);
            m.q = sol(:,2);

            % compute and set additional data
            m.p_in = p_in(m.time);
            [m.sat, m.p_cloc, tau] = computeSatP_clocTau(m, m.gamma, m.q, m.p_in);
            m.p_diff = m.p_cloc + tau .* m.q / integral(m.w, 0,1);
        end

        function m = solveDAE_explicit(m, p_in, T, gamma0, NT)
            % SOLVEDAE solves the dae model using the explicit Heun's method.
            % The solution is in the properties of the returned ModelHysteretic.
            %
            % The solution process uses Heun's method for the ode in dgamma/dt and
            % determines q using fzero.
            %
            % INPUT:
            %   p_in    - inlet pressure function (depends on t)
            %   T       - final time of simulation time [0,T] (positive)
            %   gamma0  - initial interface position (0 <= gamma0 <= 1)
            %   NT      - number of timesteps (defaault 100)

            if T <= 0.0
                error('The final time T must be positive.')
            end
            if (gamma0 < 0.0 || gamma0 > 1.0)
                error('The initial interface position gamma0 must be between 0 and 1.')
            end
            if nargin < 5
                NT = 100;
            end

            % create model
            int_fun = @(x) 3.0 ./ (m.w(x).^2 .* (m.w(x) + 3 * m.slip));
            integral01 = integral(int_fun, 0, 1);

            t = linspace(0,T, NT+1)';
            dt = T/NT;
            g = zeros(size(t));
            flux = zeros(size(t));

            g(1) = gamma0;
            tmp = integral01 + (m.M - 1) * integral(int_fun, g(1), 1);
            flux(1) = fzero(@(q) m.zeta(g(1), m.w(g(1)) * m.Ca * (p_in(t(1)) - tmp .* q)) - q ./ m.w(g(1)), 1);
            for i = 2:numel(t)
                g(i) = g(i-1) + dt * flux(i-1) / m.w(g(i-1));
                tmp = integral01 + (m.M - 1) * integral(int_fun, g(i), 1);
                flux(i) = fzero(@(q) m.zeta(g(i), m.w(g(i)) * m.Ca * (p_in(t(i)) - tmp .* q)) - q ./ m.w(g(i)), flux(i-1));
                g(i) = g(i-1) + dt * (flux(i-1) / m.w(g(i-1)) + flux(i) / m.w(g(i))) / 2;

                % check if gamma left the region [0,1]
                [value, ~, ~] = ModelHysteretic.stopAtEnds(t(i), [g(i), flux(i)]);
                if any(value < 0)
                    t(i+1:end) = [];
                    g(i+1:end) = [];
                    flux(i+1:end) = [];
                    break
                end
            end

            m.time = t;
            m.gamma = g;
            m.q = flux;

            % compute and set additional data
            m.p_in = p_in(m.time);
            [m.sat, m.p_cloc, tau] = computeSatP_clocTau(m, m.gamma, m.q, m.p_in);
            m.p_diff = m.p_cloc + tau .* m.q / integral(m.w, 0,1);
        end

        function m = solveODE(m, q, T, gamma0, options)
            % SOLVEODE solves the ode model. The solution is in the properties of the returned ModelHysteretic.
            %
            % INPUT:
            %   q       - total flux function (depends on t)
            %   T       - final time of simulation time [0,T] (positive)
            %   gamma0  - initial interface position (0 <= gamma0 <= 1)
            %   options - odeset for dae-solver options (otional)

            if T <= 0.0
                error('The final time T must be positive.')
            end
            if (gamma0 < 0.0 || gamma0 > 1.0)
                error('The initial interface position gamma0 must be between 0 and 1.')
            end
            if nargin < 5
                options = odeset();
            end
            options = odeset(options, 'Events', @ModelHysteretic.stopAtEnds);

            % solve model
            ode = @(t, y) q(t) ./ m.w(y);
            [m.time,m.gamma] = ode45(ode, [0,T], gamma0, options);

            % compute and set additional data
            m.q = q(m.time);
            [m.sat, m.p_cloc, tau] = computeSatP_clocTau(m, m.gamma, m.q);
            m.p_diff = m.p_cloc + tau .* m.q / integral(m.w, 0,1);
            int_fun = @(x) 3.0 ./ (m.w(x) .* m.w(x) .* (m.w(x) + 3 * m.slip));
            integral01 = integral(int_fun, 0, 1);
            m.p_in = m.q .* (integral01 + (m.M - 1) * arrayfun(@(g) integral(int_fun, g,1), m.gamma)) ...
                   + m.p_cloc;
        end

        function [f] = plot(m, dae)
            %PLOT plots the solution over time in one figure f.
            %
            % INPUT:
            %   dae - whether the dae model was solved (bool, ode otherwise)

            fig_title = ['solution for M = ' num2str(m.M) ', slip = ' num2str(m.slip) ', Ca = ' num2str(m.Ca)];
            if nargin < 2
                f = figure('Name', fig_title);
            elseif dae
                f = figure('Name', ['DAE ' fig_title]);
            else
                f = figure('Name', ['ODE ' fig_title]);
            end

            stackedplot(m.time, [m.gamma, m.p_in, m.q, m.sat, m.p_cloc, m.p_diff], ...
                        'DisplayLabels',["gamma", "p_in", "q", "S", "p_cloc", "p_diff"])
        end

        function [f] = plotPSat(m)
            %PLOTPSAT plots the static p_c - sat and the (solution dependent)
            % dynamic p_c - sat and p_diff - sat curves in one figure f.

            f = figure('Name', ['p_diff - sat relation for M = ' num2str(m.M) ', slip = ' num2str(m.slip) ', Ca = ' num2str(m.Ca)]);
            % dynamic p_diff - s
            plot(m.sat,m.p_diff)
            hold on
            plot(m.sat,m.p_cloc)
            % static p_c - s
            g = linspace(0,1, 100);
            % compute the saturation S = W(gamma) / W(1)
            W1 = integral(m.w, 0,1);
            s = arrayfun(@(g) integral(m.w, 0,g) ./ W1, g);
            % compute the local capillary pressure
            p_c = @(g) fzero(@(p) m.zeta(g, p), 0) ./ (m.Ca * m.w(g));
            p_c = arrayfun(p_c, g);
            plot(s, p_c)
            hold off
            legend('dyn. p_{diff}', 'dyn. p_{c,loc}', 'static p_c')
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
            fprintf(fileID,'# Ca = %g  M = %g  slip = %g  w = %s  zeta = %s\n', ...
                    m.Ca, m.M, m.slip, func2str(m.w), strrep(func2str(m.zeta), ',', ';'));
            fprintf(fileID, 'time,gamma,p_in,q,s,p_cloc,p_diff\n');
            fprintf(fileID, '%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
                            [m.time, m.gamma, m.p_in, m.q, m.sat, m.p_cloc, m.p_diff]');
            fclose(fileID);
        end

        function [] = saveCapillarityParameters(m, filename)
            %SAVECAPILLARITYPARAMETERS writes a csv table with the local static capillary pressure and the dynamic coefficient.
            %
            % INPUT:
            %   filename - the filename to use (default 'YYYY-MM-DD_HH-MM-SS.dat')

            % compute values
            g = linspace(0,1, 200)';
            [s, p_c, tau] = computeSatP_clocTau(m, g, zeros(size(g)));

            % create filename if not given
            if nargin < 2
                t = datetime();
                t.Format = 'uuuu-MM-dd_HH-mm-ss';
                filename = string(t) + '.dat';
            end

            % write csv file
            fileID = fopen(filename,'w');
            fprintf(fileID,'# Ca = %g  M = %g  slip = %g  w = %s  zeta = %s\n', ...
                    m.Ca, m.M, m.slip, func2str(m.w), strrep(func2str(m.zeta), ',', ';'));
            fprintf(fileID,'s,p_cloc,tau\n');
            fprintf(fileID,'%.6g,%.6g,%.6g\n', [s, p_c, tau]');
            fclose(fileID);
        end
    end

    methods (Access = protected)
        function [sat, p_cloc, tau] = computeSatP_clocTau(m, gamma, q, p_in)
            %COMPUTESATP_CLOCTAU computes the saturation, the local capillary pressure and the dynamic coefficient from a solution gamma, p_in, q.
            %
            % The input may be vectorial, but must have the same size.
            % If p_in is given, the local capillary pressure is computed based on p_in, instead of being estimated via zeta.
            %
            % INPUT:
            %   gamma - the interface position
            %   q     - the total flux
            %   p_in  - the inlet pressure (optional)

            % compute the saturation S = W(gamma) / W(1)
            W1 = integral(m.w, 0,1);
            sat = arrayfun(@(g) integral(m.w, 0,g) / W1, gamma);

            % compute the local capillary pressure
            if(nargin < 4)
                p_c = @(g, q) fzero(@(p) m.zeta(g, p) - q / m.w(g), 0) ./ (m.Ca * m.w(g));
                p_cloc = arrayfun(p_c, gamma, q);
            else
                int_fun = @(x) 3.0 ./ (m.w(x).^2 .* (m.w(x) + 3 * m.slip));
                integral01 = integral(int_fun, 0, 1);
                p_cloc = p_in - q .* (integral01 + (m.M - 1) * arrayfun(@(g) integral(int_fun, g,1), gamma));
            end

            % compute the dynamic coefficient
            % tau = 3 (W(1))^2 (\int_0^gamma W(x) / (w(x))^2 (3 slip + w(x))) dx / W(gamma)
            %                  + m.M * \int_gamma^1 (W(1) - W(x)) / ((w(x))^2 (3 slip + w(x))) dx / (W(1) - W(gamma)))
            tmp1 = @(x) integral(m.w, 0,x) ./ ((m.w(x)).^2 .* (3 * m.slip + m.w(x)));
            tmp2 = @(x) integral(m.w, x,1) ./ ((m.w(x)).^2 .* (3 * m.slip + m.w(x)));
            int1 = @(g) integral(@(x) arrayfun(tmp1,x), 0,g) ./ integral(m.w, 0,g);
            int2 = @(g) integral(@(x) arrayfun(tmp2,x), g,1) ./ integral(m.w, g,1);
            tau = arrayfun(@(g) 3 * W1^2 * (int1(g) + m.M * int2(g)), gamma);
            % correct for gamma == 0 or gamma == 1 (int1 == 0 or int2 == 0)
            tau(gamma == 0) = 3 * W1^2 * m.M * int2(0);
            tau(gamma == 1) = 3 * W1^2 * int1(1);
        end
    end

    methods (Static, Access = protected)
        function [value, isterminal, direction] = stopAtEnds(~,y,~)
            % STOPATENDS computes whether the simulation needs to stop early when gamma reaches 0 or 1.
            %
            % This function is used for odeXX in solveODE and solveDAEX.

            value = [y(1); 1 - y(1)];
            isterminal = ones(2,1);
            direction = -ones(2,1);
        end
    end
end
