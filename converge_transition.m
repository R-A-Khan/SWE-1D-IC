% Run data assimilation for various cases and plot results
clear all
% Parameters for the data assimilation calculations
n_obs       = [2 3 4 5 6];% Number of observation points
rand_x0     = false;      % Use random locations for observation points
x0_min      = 0.2;        % Initial observation point
ntrial      = 1;          % Number of trials
N           = 1024;       % Number of grid points
tmax        = 2;          % control time
xmax        = 3;          % domain size: xmin = -xmax
iter_max    = 1000;        % Max number of iterations
line_min    = true;       % Find optimal step
smooth_grad = false;      % Smooth L2 gradient to Sobolev gradient
filt        = 0.01;        % filtering parameter for Sobolev gradient   
% conj_grad_type = 0;       % 0 = Steepest Descent, 1 = Fletcher Reeves, 2 = Polak-Ribière
% iter_chunk = 5;           % No. of iterations after which descent algorithm resets to steepest descent 
% noisy_obs  = false;       % True = noisy observations and includes tikhonov regularisation
% tikh_par   = 0.01;        % Tikhonov Regularisation Parameter
% reg_norm   =  1;          % Type of norm for regularisation: L2 =1, H1 (seminorm = 22) = 21, H2 = 31 (seminorm = 32)
bath       = true;        % Include non-flat bottom bathymetry
bath_func = @(x) 0.05.*exp(-(x-1.5).^2/0.1^2); % function for bottom topography
%     bath_func = @(x) 0.05.*(tanh(2*(x+2)) - tanh(2*(x-2)) );

% dmu = 0:0.0001:0.25; % Spacing between observation points
% dmu = 2*xmax/N*4:0.01:0.5; % Spacing between observation points

% L=2*pi; dk = 2*pi/6;
% %% Condition on psi(k)
% %dk = 1;
% dk = 0.075;
% k=dk*(0:1:N/2);
% kmax = 30; % Maximum wavenumber of initial conditions
% 
% 
% for jk = 1:numel(k)
%     for j_obs = n_obs
%         j = find(n_obs==j_obs);
%         for l = 1:numel(dmu)
%             x0_max = x0_min+(j_obs-1)*dmu(l);
%             x0_vals = x0_min:dmu(l):x0_max;
%             error(l,jk,j) = abs(1/j_obs * sum(exp(1i*2*k(jk)*x0_vals)));
%         end
%     end
%     plot(dmu,error(:,jk,1),'k-');hold on;
%     if k(jk)>kmax
%         break
%     end
%     %plot(dmu,err(:,jk,2),'r-');hold on;
% end
% %plot(dmu,err(:,1),'k-');hold on;
% %plot(dmu,err(:,2),'k--');
% %plot(dmu,err(:,3),'k-.');
% %legend('N_{obs}=2','N_{obs}=4','N_{obs}=6');
% xlabel('$\Delta x$','Interpreter','Latex','fontsize',20);
% ylabel('$|\hat{\psi}(k)|$','Interpreter','Latex','fontsize',20);grid on;
% axis([0.1 0.25 0.99 1]);
% set(gca,'FontSize',20); grid on; 
% hold off;
%%
% print -depsc2 ~/doc/papers/data_assimilation/figures/N_2_k_0_30_dk_0075.eps

%%
dmu = [0.09 0.375 0.48];
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
        for l = 1:numel(dmu)
            disp(['spacing = ', sprintf('%0.3e',dmu(l))]);
            [eta_optimum(:,k,l,j_nl), eta_exct0(:,k,l,j_nl),bathymetry(:,k,l,j_nl), X, err(:,k,l,j_nl), err_std(:,k,l,j_nl),...
            grad(:,k,l,j_nl), grad_std(:,k,l,j_nl), cost(:,k,l,j_nl), cost_std, n_iter(k,l,j_nl)] = ...
            data_assimil(nl,N,j_obs,rand_x0,x0_min,dmu,ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, bath, bath_func);
        end
    end
    disp(' ');
end
% Save results
save('IC_Assimilation_with_bathy', 'ntrial', 'n_obs', 'rand_x0','x0_min', 'dmu', 'line_min','smooth_grad', 'iter_max', ...
    'eta_optimum', 'bathymetry', 'eta_exct0', 'X', 'err', 'err_std', 'grad', 'grad_std', 'cost', 'cost_std', 'n_iter');
% %% Plotting
% clear all;
% %load('~/doc/papers/data_assimilation/data/Assimilation_data_transition_cost');
% load('Assimilation_data_transition');
% clf;
% iter = 0:iter_max;
% marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
% marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};
% 
% skip=iter_max/25;
% 
% n_dmu = numel(dmu);
% 
% figure(2);
% for j_obs = [2]
%     k = find(n_obs==j_obs)
%     for j_nl = 1:2
%         n(:,k,j_nl) = n_iter(k,:,j_nl);
%     end
%     for l = 1:n_dmu
%         for j_nl = 1:2
%             err_final(l,j_nl) = err(n(l,k,j_nl),k,l,j_nl);
%         end
%     end
%     h(k)=plot(dmu, err_final(:,1),'k-','linewidth',1.5);hold on;
%     plot(dmu, err_final(:,2),'k--','markersize',10,'linewidth',1.5);
%     for l = 1:n_dmu
%         if n_iter(k,l,1)<5001
%             plot(dmu(l),err_final(l,1),'k','markersize',10);
%         else
%             plot(dmu(l),err_final(l,1),'k','markersize',10,'MarkerFaceColor','k');
%         end
%         if n_iter(k,l,2)<5001
%             plot(dmu(l),err_final(l,2),'k--','markersize',10);
%         else
%             plot(dmu(l),err_final(l,2),'k--','markersize',10,'MarkerFaceColor','k');
%         end
%     end
% end
% set(gca,'FontSize',16); grid on; 
% %legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
% %legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4'});
% legend('Linear','Nonlinear');
% xlabel('$\Delta x$','interpreter','latex','fontsize',24); 
% % ylabel('$||f-\tilde{\phi}||_2/||f||_2$','interpreter','latex','fontsize',24);xlim([0 0.5]);ylim([0 0.35]);
% % print -depsc2 ~/doc/papers/data_assimilation/figures/converge_transition.eps
% hold off;
% %%
% plot(X,eta_optimum(:,1,1,1),'k-.','linewidth',1.5);grid on; hold on;
% plot(X,eta_optimum(:,1,2,2),'k--','linewidth',1.5);
% plot(X,eta_optimum(:,1,2,2),'k-','linewidth',1.5);
% 
% set(gca,'FontSize',16);
% legend('Linear','Nonlinear');
% xlabel('x'); ylabel('\eta_0(x)');xlim([min(X) max(X)]);ylim([-0.05 0.05])
% % print -depsc2 optimum_eta.eps
% hold off; 