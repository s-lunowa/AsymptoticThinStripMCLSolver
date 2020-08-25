function example()
% EXAMPLE shows how to solve the DAE / ODE of the Model with dynamic contact angle model.
% The solutions are plotted and saved for further usage.
%
% (c) 2020 Stephan B. Lunowa
%
% This work is licensed under the Creative Commons Attribution 4.0 International License.
% You should have obtained a LICENCE file alongside this file.
% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

% The model parameters need to be given first:

% The wall function defines the shape of the thin strip.
% It must be positive, vectorized and satisfy w(0) = 1.
% Here we use a 'constriction'
w = @(x) 2.0/3.0 + cos(2 * pi * x) / 3.0;

% The contact angle model is given in terms of the position and the speed
% of the contact line. It must attain values in [0,pi], and must be vectorized.
% Here we use a linear structure in the position and the velocity, restricted
% to [0,pi]. Note that this has no physical meaning.
theta = @(x, u) max(min((x + 1 - u) * (pi / 6), pi), 0);
% The contact angle can also be static, e.g. theta = @(x,u) ones(size(u))

% The effective capillary number must be positive.
Ca = 0.5;

% The viscosity ratio of fluid II to fluid I must be non-negative.
% A value of zero results in one-phase dynamics for fluid I.
M = 0.5;

% The slip length at the wall must be non-negative.
% A value of zero indicates only slip at the interface, while values above
% zero result in slip everywhere.
slip = 0;

% With these parameters, we can create the model.
m = Model(w, theta, Ca, M, slip);

% At this point, we can already compute the static capillaty pressure - saturation
% and dynamic coefficient - saturation curves and save them into the file 'example_pc-tau-sat.dat'.
m.saveCapillarityParameters('example_pc-tau-sat.dat')

% To simulate the solution, we need to set a positive final time.
T = 1;

% We also need an initial position of the interface, which must be within [0,1].
% Note that the simulation will stop immediately if gamma reaches 0 or 1.
gamma0 = 1e-3;

% The ODE model needs a given total flux function.
% Here we choose a constant one for simplicity.
q = @(t) ones(size(t));

% The DAE model needs a given inlet pressure function.
% Here we choose a constant one for simplicity.
p_in = @(t) 12 * ones(size(t));

% Finally, we can set typical options for the solvers such as the absolute
% and relative tolerance, or the maximal time-step size.
% Note that the DAE solver has two components, so that it needs two
% absolute tolerances, while the ODE solver only needs one.
optionsDAE = odeset('MaxStep', 1e-2, 'RelTol',1e-3, 'AbsTol',[1e-6, 1e-7]);
optionsODE = odeset(optionsDAE, 'AbsTol',1e-6);

% Solving the models is done by the solveDAE / solveODE functions.
mDAE = m.solveDAE(p_in, T, gamma0, optionsDAE);
mODE = m.solveODE(q, T, gamma0, optionsODE);

% We save the solutions into the file 'example_solutionXXX.dat'.
mDAE.saveSolution('example_solutionDAE.dat')
mODE.saveSolution('example_solutionODE.dat')

% The solutions can be easily plotted with a lot of additional information.
% The bool argument only determines if DAE or ODE appears in the title.
f1 = mDAE.plot(true);
f2 = mODE.plot(false);

% you need to manually continue here
pause
if(isvalid(f1)); close(f1); end
if(isvalid(f2)); close(f2); end

% The pressure curves of the solutions can also be plotted against the saturation.
f1 = mDAE.plotPSat();
f2 = mODE.plotPSat();

% you need to manually continue here
pause
if(isvalid(f1)); close(f1); end
if(isvalid(f2)); close(f2); end

end