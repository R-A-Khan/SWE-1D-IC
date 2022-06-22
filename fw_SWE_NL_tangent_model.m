function trend = fw_SWE_NL_tangent_model(t, y, h, u, eta)
% Usage: trend = fw_SWE_NL_periodic_v2(t, y, h)
%
% Finds trend of nonlinear SWE equations with periodic boundary conditions
% using second order finite difference / finite volume approximation
% grid size is 2*N, with dh/dt evaluated at centres and du/dt at edges
%
% 1 to N rows of dy represent     dh/dt = - ((1+eta)*u)_x
% N+1 to 2*N rows of dy represent du/dt = - (0.5*u^2 + eta)_x 
%
% Input:
% t = time (not used)
% y = input variables (height and velocity in one vector)
% h = spacing of grid points in space
%
% Output:
% trend = trend for both height and velocity equations in one vector

N = length(y)/2;
eta_hat = y(1:N);
u_hat   = y(N+1:2*N);

trend = zeros(2*N, 1);
trend(1:N)     =  -( cent_diff_h(h, u.*mid2edge(eta_hat)) + cent_diff_h(h, (1+mid2edge(eta)).*u_hat) );
trend(N+1:2*N) =  -( cent_diff_u(h,edge2mid(u.*u_hat)) + cent_diff_u(h,eta_hat) );
