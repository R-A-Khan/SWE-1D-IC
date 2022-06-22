function J = line_min_mult(obs_x0_mult, eta0_adj, N, h, x0_inds, tau, eta0, dt, tmin, tmax, funsolve, beta)
% Usage:  = line_min_mult(obs_x0_mult, eta0_adj, N, h, x0_inds, tau, eta0, dt, tmin, tmax, funsolve)
%
% Calculates cost function for line minimization algorithm
%
% Input:
% obs_eta_t = observation vector at x0 for all t
% eta0_adj  = adjoint height at t0 for all x
% tau = stepsize for line minimisation
% eta0 = initial condition at t0 for all x
% dt = time step size
% [tmin, tmax] = time range
% funsolve = function for FD spatial discretisation of forward solver
%
% Output:
% J = cost function

N0 = eta0 + tau.*eta0_adj;

u = zeros(N,1);
H = [N0;u];
[~, ~, eta_x0_t_mult, ~, T1, ~] = FW_solve(h, dt, H, tmin, tmax, x0_inds, funsolve, beta);

% Define Integrand
sum = zeros(1,length(T1));
for i = 1:length(x0_inds)
    sum = sum + (obs_x0_mult(i,:) - eta_x0_t_mult(i,:)).^2 ;
end
% Numerical Integration using trapezoid Rule
J = trapz(T1,0.5*sum);

% tau_n is argmin of J over tau
% value of tau  where J is a minimum.
end

