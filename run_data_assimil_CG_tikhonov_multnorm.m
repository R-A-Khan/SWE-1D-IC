% Run data assimilation for various cases with Tikhonoc regularisation and plot results
clear all
% Parameters for the data assimilation calculations
n_obs       = [5];        % Number of observation  points
rand_x0     = false;      % Use random locations for observation points
x0_min      = 0.2;        % Location of first observation point 
dmu         = 0.2;        % Spacing between observation points
ntrial      = 1;          % Number of trials
N           = 1024;       % Number of grid points
tmax        = 2;          % control time 
xmax        = 3;          % domain size: xmin = -xmax
iter_max    = 500;        % Max number of iterations
line_min    = true;       % Find optimal step
smooth_grad = false;      % Smooth L2 gradient to Sobolev gradient
filt        = 0.1;        % filtering parameter for Sobolev gradient   
conj_grad_type = 2;       % 1 = Fletcher Reeves, 2 = Polak-Ribière
iter_chunk = 5;           % No. of iterations after which descent algorithm resets to steepest descent 
noisy_obs  = false;       % True = noisy observations and includes tikhonov regularisation
tikh_par   = 0.01;        % Tikhonov Regularisation Parameter
reg_norm   =  1;          % Type of norm: L2 =1, H1 (H1 seminorm = 22) = 21, H2 = 31 ((H2 seminorm = 32)
bath       = true;        % Include non-flat bottom bathymetry


disp(' ');
for j_nl=1:2
    if (j_nl==1)
        disp('Linear SWE');
    else
        disp('Nonlinear SWE');
    end
    nl=j_nl-1;
    for j_obs = n_obs
        k = find(n_obs==j_obs);
%         dmu  = abs((xmax-1e-4) - x0_min)/n_obs(k);
        [eta_optimum(:,k,j_nl), eta_exct0(:,k,j_nl),bathymetry(:,k,j_nl), X, err(:,k,j_nl), err_std(:,k,j_nl), grad(:,k,j_nl), grad_std(:,k,j_nl), cost(:,k,j_nl), cost_std] = ...
            data_assimil_CG_tikhonov_multnorm(nl,N,j_obs,rand_x0,x0_min,dmu,ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, conj_grad_type, iter_chunk, noisy_obs, tikh_par, reg_norm, bath);
    end
    disp(' ');
end
% Save results
save('Assimilation', 'ntrial', 'n_obs', 'rand_x0','x0_min', 'dmu', 'line_min','smooth_grad', 'iter_max', ...
    'eta_optimum', 'bathymetry', 'eta_exct0', 'X', 'err', 'err_std', 'grad', 'grad_std', 'cost', 'cost_std');
%% Plotting
clear all;
%

load('Assimilation');
clf;
iter = 0:iter_max;
marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/25;

set(0,'DefaultFigureWindowStyle','docked')

figure(1);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k)=semilogy(iter(1:skip:end), err(1:skip:end,k,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), err(1:skip:end,k,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
% legend(h, {'N_{obs} = 2','N_{obs} = 3', 'N_{obs} = 4'});
xlabel('Iteration n'); ylabel('$|\eta_0^{(n)}(x) - f(x)|_2/|f(x)|_2$','interpreter','latex');xlim([0 iter_max]);
ylim([1e-6 1e0]);
print -depsc2 err_assimil.eps
hold off;

figure(2);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k)=semilogy(iter(1:skip:end), grad(1:skip:end,k,1)/grad(1,k,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), grad(1:skip:end,k,2)/grad(1,k,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$|\nabla J^{(n)}(x)|_2/|\nabla J^{(0)}(x)|_2$','interpreter','latex');xlim([0 iter_max]);
print -depsc2 err_grad.eps
hold off;

figure(3);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k) = semilogy(iter(1:skip:end), cost(1:skip:end,k,1)/cost(1,k,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
           semilogy(iter(1:skip:end), cost(1:skip:end,k,2)/cost(1,k,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$J^{(n)}/J^{(0)}$','interpreter','latex');xlim([0 iter_max]);
print -depsc2 err_cost.eps
hold off;

figure(4);
plot(X,eta_optimum(:,end,1),'k-','linewidth',1.5);grid on; hold on;
plot(X,eta_optimum(:,end,2),'k--','linewidth',1.5);
set(gca,'FontSize',16);
legend('Linear','Nonlinear');
xlabel('x'); ylabel('\eta_0(x)');xlim([min(X) max(X)]);
print -depsc2 optimum_eta.eps
hold off;

figure(5); plot(X,eta_exct0(:,end,1), 'linewidth',2); hold on
plot(X, (bathymetry(:,end,1)-1), 'linewidth', 2);hold off

% figure(5)
% for j_obs = n_obs
%     k = find(n_obs==j_obs);
%     semilogy(j_obs, err_std(end,k,1)/err(end,k,1)*100, marker1{k},'markersize',10,'linewidth',1.5);hold on;
%     semilogy(j_obs, err_std(end,k,2)/err(end,k,2)*100, marker2{k},'markersize',10,'linewidth',1.5);
% end
% xlabel('N_{obs}');ylabel('% error');
% for j_obs = n_obs
%     k = find(n_obs==j_obs);
%     disp([sprintf('%1.f',j_obs),' ',sprintf('%0.3e',err(end,k,1)), ' ', sprintf('%0.3e',err_std(k,1)), ...
%        ' ',sprintf('%0.3e',err(end,k,2)), ' ', sprintf('%0.3e',err_std(k,2))]);
% end

% figure(5); 
% semilogy(X, abs(eta_optimum(:,1)-eta_exct0(:,1)), 'k-');hold on;
% semilogy(X, abs(eta_optimum(:,2)-eta_exct0(:,2)), 'k--');
% set(gca,'FontSize',16); grid on;
% xlabel('x'); ylabel('Error'); xlim([-xmax xmax]);
