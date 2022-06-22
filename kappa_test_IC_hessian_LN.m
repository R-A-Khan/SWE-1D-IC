clear all
% Calculates kappa test for linear or nonlinear SWE for a range of
% perturbation sizes

% Parameters for the data assimilation calculations
nl          = false;      % nonlinear or linear equations
n_obs       = 3;          % Number of observation points
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



    fw = @fw_SWE_periodic;
    bw = @bw_SWE_periodic_mult;
    fw_tangent = @fw_SWE_NL_tangent_model_LN;
    bw_SOA = @bw_SWE_SOA_LN;

    
    
    % Exact initial conditions
    amp = 0.05; eta_exct0 = amp.*exp(-X.^2/0.1^2);
    [x0_inds, x0_pts] = mult_x0(X,x0_vals);
    
    % STEP I: Run forward solver using exact initial conditions to get observations
    u = zeros(N,1);
    H = [eta_exct0 ; u];
    [~, ~, obs_x0_mult, obs_T_x, T_hat, ~] = FW_solve(h, sz, H, 0, tmax, x0_inds, fw, beta);
    
    
    % STEP 2: Distort exact IC to get "guess" IC eta0(:,1) used in algorithm with analytical grad_J
    eta0 = 0.9*eta_exct0;
    % Perturbation direction
    eta0_prime_g = amp*cos(X);
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
%   Hessian_J = zeros(N,N);
 
    % Solve tangent model for perturbed system with eta_hat(x,0) = e_k 
%  for k = 1:N
        u = zeros(N,1);
        eta_hat0 = amp*sin(X);
%         eta_hat0(k) = 1;
        H_hat = [eta_hat0 ; u];
        [~, ~, eta_hat_x0_pts, test_l, T_hat, Y_hat] = FW_solve_tangent_LN(h, sz, H_hat, 0, tmax, x0_inds, fw_tangent, beta);
        eta_hat = Y_hat(1:N,:);
        u_hat = Y_hat(N+1:2*N,:);

        % Solve SOA model for eta_bat(x,0)
        H_bar = zeros(2*N,1);
    %     [eta_bar_t0_x, T_bar, Y_bar] = BW_solve_SOA_LN(h, sz, H_bar, 0, tmax, x0_inds, bw_SOA, eta_hat_x0_pts);
        [eta_bar_t0_x, T_bar, Y_bar] = BW_solve(h, sz, H_bar, 0, tmax, x0_inds, bw_SOA, zeros(length(x0_inds),length(T_hat)), eta_hat_x0_pts, u_all, eta_all, T_hat, T2, nl);
        % Define Hessian

        Hessian_J =  - eta_bar_t0_x;
%  end
%     plot(-eta_bar_t0_x)
% % %     
%     figure(6);
%     for i = 1:length(T_hat)
%         % Example of plot
%         Z = eta_hat(:,i);
%         plot(X,Z, 'r');
%         ylim([-1e-1 1e-1])
%         xlim([xmin xmax])
%         M(i)=getframe(gcf);
% 
%     end


tau = 10^(-1); 
ep = logspace(-8,4,200);
for i = 1:length(ep)
    % STEP 1(iv)
    % Compute eta_pert at the point x0 given perturbed IC eta + ep*eta_prime
    eta_pert = eta0 + ep(i)*eta0_prime_g;
    
    eta_pert_tau = eta0 + tau*eta0_prime_g ;
    eta_pert_ep = eta0 + ep(i)*eta_hat0 ;
    eta_pert_ep_tau = eta0 + ep(i)*eta_hat0 + tau*eta0_prime_g;

    % Run forward solver to get eta_t at observation points given "initial guess" IC eta_0
    up = zeros(N,1);
    Hp = [eta_pert;up];
    [~, ~, etapert_x0_t_mult, etapert_T_x, Tp, Yp] = FW_solve(h, sz, Hp, 0, tmax, x0_inds, fw, beta);
    
    H_tau = [ eta_pert_tau; up];
    H_ep =  [eta_pert_ep; up];
    H_tau_ep = [eta_pert_ep_tau; up];
    
    [~, ~, etapert_x0_t_tau, ~, ~,~] = FW_solve(h, sz, H_tau, 0, tmax, x0_inds, fw, beta);
    [~, ~, etapert_x0_t_ep, ~, ~,~] = FW_solve(h, sz, H_ep, 0, tmax, x0_inds, fw, beta);
    [~, ~, etapert_x0_t_tau_ep, ~, ~,~] = FW_solve(h, sz, H_tau_ep, 0, tmax, x0_inds, fw, beta);
    
    
    
    % (J'(phi + ep*eta''; \eta')
    sum_1 = zeros(1,length(T2));
    for j = 1:length(x0_inds)
        sum_1 = sum_1 + (obs_x0_mult(j,:) - etapert_x0_t_tau_ep(j,:)).^2 ;
    end
    J_1 = trapz(T2, 0.5*sum_1);
    
   sum_2 = zeros(1,length(T2));
    for j = 1:length(x0_inds)
        sum_2 = sum_2 + (obs_x0_mult(j,:) - etapert_x0_t_ep(j,:)).^2 ;
    end
    J_2 = trapz(T2, 0.5*sum_2);
    
   J_prime_1 = (1/tau)*(J_1-J_2); 
    

    % (J'(phi'; \eta')
    sum_3 = zeros(1,length(T2));
    for j = 1:length(x0_inds)
        sum_3 = sum_3 + (obs_x0_mult(j,:) - etapert_x0_t_tau(j,:)).^2 ;
    end
    J_3 = trapz(T2, 0.5*sum_3);
    
   sum_4 = zeros(1,length(T2));
    for j = 1:length(x0_inds)
        sum_4 = sum_4 + (obs_x0_mult(j,:) - eta_x0_t_mult(j,:)).^2 ;
    end
    J_4 = trapz(T2, 0.5*sum_4);
    
   J_prime_2 = (1/tau)*(J_3-J_4); 
   
   % Trying \int -l^l - eta* eta'dx definition of J_prime_1 and J_prime_2
   [eta_a_pert_t0_x, T3, Y3] = BW_solve(h, sz, H3, 0, tmax, x0_inds, bw, obs_x0_mult, etapert_x0_t_ep, u_all, eta_all, T_hat, T2, nl);
    
   Riesz_grad_1 = trapz(X, (-eta_a_pert_t0_x).*eta0_prime_g);

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

    Gateaux_grad = (1/ep(i))*(J_eta_pert-J_eta);
    Riesz_grad = trapz(X, grad_J.*eta0_prime_g);
    
    

    kappa(i) = Gateaux_grad/Riesz_grad;

    psi_equiv_1 = (1/ep(i))*(J_prime_1 - J_prime_2);
%     psi_equiv_2 = (0.5*ep(i)^2)* trapz(X,(trapz(X,eta0_prime_g.*Hessian_J)) .* eta0_prime_g)
    psi_equiv_2 =  trapz(X, (- eta_bar_t0_x.* eta0_prime_g ));
    
    psi_equiv_1b = (1/ep(i))*(Riesz_grad_1 - Riesz_grad);
    
      
    psi(i) = psi_equiv_1/psi_equiv_2;
    psi_2(i) = psi_equiv_1b/psi_equiv_2;
end
    
% figure(1);
% semilogx((ep(1:length(kappa))),abs(kappa),'k',  'linewidth',1.5); grid on;
% set(gca,'FontSize',22)
% legend('Nonlinear SWE', 'location','northeastoutside');
% %title('Kappa test for Numerical Solution I')
% xlabel('log(\epsilon)')
% ylabel('|\kappa(\epsilon)|');
% %%
% 
% 
% %%
% figure(2);
% loglog(ep(1:length(kappa)), abs(abs(kappa)-1),'k',  'linewidth',1.5); grid on;
% set(gca,'FontSize',22)
% %title('Kappa test for Numerical  Solution  II')
% legend('Nonlinear SWE', 'location','northeastoutside');
% xlabel('log(\epsilon)');
% ylabel('log|\kappa(\epsilon)-1|');
%%

figure(3);
semilogx((ep(1:length(psi))),abs(psi),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
legend('Nonlinear SWE', 'location','northeastoutside');
%title('Kappa test for Numerical Solution I')
xlabel('log(\epsilon)')
ylabel('|\psi(\epsilon)|');

figure(4);
loglog(ep(1:length(psi)), abs(abs(psi)-1),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
%title('Kappa test for Numerical  Solution  II')
legend('Nonlinear SWE', 'location','northeastoutside');
xlabel('log(\epsilon)');
ylabel('log|\kappa(\epsilon)-1|');

figure(5);
semilogx((ep(1:length(psi_2))),abs(psi_2),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
legend('Nonlinear SWE', 'location','northeastoutside');
%title('Kappa test for Numerical Solution I')
xlabel('log(\epsilon)')
ylabel('|\psi_b(\epsilon)|');
%%


figure(6);
loglog(ep(1:length(psi_2)), abs(abs(psi_2)-1),'k',  'linewidth',1.5); grid on;
set(gca,'FontSize',22)
%title('Kappa test for Numerical  Solution  II')
legend('Nonlinear SWE', 'location','northeastoutside');
xlabel('log(\epsilon)');
ylabel('log|\psi_b(\epsilon)-1|');