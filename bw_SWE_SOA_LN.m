function trend = bw_SWE_SOA_LN(t, y, h, x0_inds, obs, eta_hat_x0, ~, ~)
% Usage: trend = bw_SWE_periodic_mult(t, y, h, x0_inds, obs, eta)
%
% Finds trend of backwards (adjoint) linear SWE
% using second order finite difference / finite volume approximation
% grid size is 2*N, with dh/dt evaluated at centres and du/dt at edges
%
% Input:
% t       = time (not used)
% y       = input variables (height and velocity in one vector)
% N       = number of grid points in space
% h       = spacing of grid points
% x0_inds = indices of observation positions
% obs     = height at observation position from t=tmin to t=tmax
% eta     = height from forward solver using approximate initial conditions
%
% Output:
% trend = trend for both height and velocity equations in one vector

N = length(y)/2;
eta_bar = y(1:N);
u_bar   = y(N+1:2*N);

trend = zeros(2*N, 1);
trend(1:N)     = cent_diff_h(h, u_bar); 
trend(N+1:2*N) = cent_diff_u(h, eta_bar);

% Add observations
for i = x0_inds
    k = find(x0_inds==i);
    trend(i) = trend(i) - eta_hat_x0(k)/h;
end
