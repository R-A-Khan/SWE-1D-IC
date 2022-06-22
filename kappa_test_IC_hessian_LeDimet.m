
clear all
% Calculates kappa test for linear or nonlinear SWE for a range of
% perturbation sizes
 
% Parameters for the data assimilation calculations
nl          = true;      % nonlinear or linear equations
n_obs       = 6;          % Number of observation points
x0_min      = 0.2;        % Location of first observation point
dmu         = 0.2;       % Spacing between observation points
N           = 256;        % Number of grid points
tmax        = 2;          % control time
xmax        = 3;          % domain size: xmin = -xmax
smooth_grad = false;      % Smooth L2 gradient to Sobolev gradient
filt        = 0.4;       % filtering parameter for Sobolev gradient
bath        = false;
% Observation points
x0_max = x0_min+(n_obs-1)*dmu;
x0_vals = x0_min:dmu:x0_max;
 
xmin = -xmax;
h = (xmax - xmin)/N;
sz = h;
tmin = 0;
 
numSteps = round(abs((tmax-tmin)/sz));
tmax = numSteps*sz;
X = xmin:h:xmax-h;
X = reshape(X, [N,1]);
 
if bath
%     beta = b_amp.*(tanh(2*(X+2*ones(N,1))) - tanh(2*(X-2*ones(N,1))) );
      beta = bath_func(X);
else
    beta = zeros(N,1);
end
 
 
 
    fw = @fw_SWE_NL_periodic_with_bath;
    bw = @bw_SWE_periodic_mult_NL_v2;
    fw_tangent = @fw_SWE_NL_tangent_model;
    bw_SOA = @bw_SWE_SOA;
 
     
    % STEP 1(iii)
    % Set Initial conditions
     
    % Adding noise
    %     for i = 1:N
    %         if eta_exct0(i) > 0
    %         eta0_noisy(i) = awgn(eta_exct0(i),20,'measured');
    %         else eta0_noisy(i) = eta_exct0(i);
    %         end
    %     end
    %     save('noisy_IC_NL.mat','eta0_noisy');
     
    % Load presaved distorted IC file from above
    % load('noisy_IC_NL.mat');
    %
    % filt_h = 0.4;
    % eta0_smooth = ffft_smoothing(eta0_noisy,filt_h);
    % eta0 = eta0_smooth;
     
    % Exact initial conditions
    amp = 0.05; eta_exct0 = amp.*exp(-X.^2/0.1^2);
    [x0_inds, x0_pts] = mult_x0(X,x0_vals);
     
    % STEP I: Run forward solver using exact initial conditions to get observations
    u = zeros(N,1);
    H = [eta_exct0 ; u];
    [~, ~, obs_x0_mult, obs_T_x, T_hat, ~] = FW_solve(h, sz, H, 0, tmax, x0_inds, fw, beta);
     
     
    % STEP 2: Distort exact IC to get "guess" IC eta0(:,1) used in algorithm with analytical grad_J
    eta0 = amp*sin(X);
    % Perturbation direction
    eta0_prime_g = eta0;
    eta_prime_x0 = eta0_prime_g(x0_inds);
     
    % STEP 2.2: Run forward and backwards solvers to get eta_a_t0_x
    u2 = zeros(N,1);
    H2 = [eta0; u2];
    [eta_all, u_all, eta_x0_t_mult, eta_T_x, T2, Y2] = FW_solve(h, sz, H2, 0, tmax, x0_inds, fw, beta);
     
    H3 = zeros(2*N,1);
    [eta_a_t0_x, T3, Y3] = BW_solve(h, sz, H3, 0, tmax, x0_inds, bw, obs_x0_mult, eta_x0_t_mult, u_all, eta_all, T_hat, T2, nl);
    eta_adj = Y3(1:N,:);
    u_adj = Y3(N+1:2*N,:);
     
    % Define grad_J
    grad_J = -eta_a_t0_x;
     
     
  %Pre-Allocate Hessian 
  Hessian_J = zeros(N,N);
  for k = 1:N 
    % Solve tangent model for perturbed system with eta_hat(x,0) = e_k 
     
    u = zeros(N,1);
    eta_hat0 = zeros(N,1);
    eta_hat0(k) = 1;
    H_hat = [eta_hat0 ; u];
    [~, ~, eta_hat_x0_pts, ~, T_hat, Y_hat] = FW_solve_tangent(h, sz, H_hat, 0, tmax, u_all, eta_all, x0_inds, fw_tangent);
    eta_hat = Y_hat(1:N,:);
    u_hat = Y_hat(N+1:2*N,:);
     
    % Solve SOA model for eta_bat(x,0)
    H_bar = zeros(2*N,1);
    [eta_bar_t0_x, T_bar, Y_bar] = BW_solve_SOA(h, sz, H_bar, 0, tmax, x0_inds, bw_SOA, eta_hat_x0_pts, u_all, eta_all, u_adj, eta_adj, u_hat, eta_hat, T_hat, T2, nl);
     
    % Define Hessian
     
    Hessian_J(:,k) =  - eta_bar_t0_x;
  end
     
    % Smooth gradient
    if smooth_grad
        grad_J = grad_smooth(grad_J,filt,X);
    end
  
    ep = logspace(-8,-2,100);
    for i = 1:length(ep)
        % STEP 1(iv)
        % Compute eta_pert at the point x0 given perturbed IC eta + ep*eta_prime
        eta_pert = eta0 + ep(i)*eta0_prime_g;
         
        % Run forward solver to get eta_t at observation points given "initial guess" IC eta_0
        up = zeros(N,1);
        Hp = [eta_pert;up];
        [~, ~, etapert_x0_t_mult, etapert_T_x, Tp, Yp] = FW_solve(h, sz, Hp, 0, tmax, x0_inds, fw, beta);
         
        % STEP 2: Compute Gateaux derivative representation of gradient
        sum = zeros(1,length(T2));
        for j = 1:length(x0_inds)
            sum = sum + (obs_x0_mult(j,:) - eta_x0_t_mult(j,:)).^2 ;
        end
        J_eta = trapz(T2, 0.5*sum);
         
        sum_p = zeros(1,length(Tp));
        for j = 1:length(x0_inds)
            sum_p = sum_p + (obs_x0_mult(j,:) - etapert_x0_t_mult(j,:)).^2 ;
        end
        J_eta_pert = (trapz(Tp,0.5*sum_p));
         
        Gateaux_grad = (1/ep(i))*(trapz(Tp,0.5*sum_p)-J_eta);
        Riesz_grad = trapz(X, grad_J.*eta0_prime_g);
        scalar_prod_grad = ((trapz((obs_x0_mult - eta_x0_t_mult), 2).*eta_prime_x0));
         
         
        kappa_equiv_1 = J_eta_pert-J_eta;
        kappa_equiv_2 = ep(i).*Riesz_grad;
         
        psi_equiv_1 = J_eta_pert-J_eta - ep(i).*Riesz_grad
        psi_equiv_2 = (0.5*ep(i)^2)* trapz(X,((eta0_prime_g)'*Hessian_J)' .* eta0_prime_g)
        
         
        % STEP 4: Kappa test for analytical solution
        psi(i) = psi_equiv_1/psi_equiv_2;
        kappa(i) = Gateaux_grad/Riesz_grad;
 
         
    end
%%
figure(1);
semilogx((ep(1:length(kappa))),abs(kappa),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
legend('Nonlinear SWE', 'location','northeastoutside');
%title('Kappa test for Numerical Solution I')
xlabel('log(\epsilon)')
ylabel('|\kappa(\epsilon)|');
%%
 
 
%%
figure(2);
loglog(ep(1:length(kappa)), abs(abs(kappa)-1),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
%title('Kappa test for Numerical  Solution  II')
legend('Nonlinear SWE', 'location','northeastoutside');
xlabel('log(\epsilon)');
ylabel('log|\kappa(\epsilon)-1|');
%%
 
figure(3);
semilogx((ep(1:length(psi))),abs(psi),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
legend('Nonlinear SWE', 'location','northeastoutside');
%title('Kappa test for Numerical Solution I')
xlabel('log(\epsilon)')
ylabel('|\psi(\epsilon)|');