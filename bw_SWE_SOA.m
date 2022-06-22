function trend = bw_SWE_SOA(t, y, h, x0_inds, eta_hat_x0, u_all, eta_all,u_adj, eta_adj, u_hat, eta_hat)
% Usage: trend = bw_SWE_periodic_mult_NL_v2(t, y, h, x0_inds, obs, eta, u_all, eta_all)
%
% Finds trend of backwards (adjoint) nonlinear SWE
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
% u_all   = velocity at all positions and times
% eta_all = height at all positions and times
%
% Output:
% trend = trend for both height and velocity equations in one vector

N = length(y)/2;
eta_bar = y(1:N);
u_bar   = y(N+1:2*N);

trend = zeros(2*N, 1);
trend(1:N)     = - (edge2mid(u_all).*cent_diff_h(h,mid2edge(eta_bar)) + edge2mid(u_hat).*cent_diff_h(h, mid2edge(eta_adj)) + cent_diff_h(h, u_bar) );
trend(N+1:2*N) = - (mid2edge(eta_hat).*cent_diff_u(h,eta_adj)+((mid2edge(eta_all)+1).*cent_diff_u(h,eta_bar)) + u_all.*cent_diff_u(h, edge2mid(u_bar)) + (u_hat).*cent_diff_u(h,edge2mid(u_adj)) );

% Add observations
for i = x0_inds
    k = find(x0_inds==i);
    trend(i) = trend(i)  - eta_hat_x0(k)/h;
end


