
load('IC_Assimilation_with_bathy_amp_0_05.mat');
clf;
iter = 0:iter_max;
% marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker1 = {'k-', 'b-', 'g-', 'r-', 'y-', 'c-d', 'm-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
% marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};
marker2 = {'k*','b*','g*','r*','y*','c*','m*','kv','k>','k<','kp','kh'};


skip=iter_max/25;


obs_pts = 4;
m = find(n_obs==obs_pts);

delta_x = 0.375;
l =find(dmu==delta_x);



figure(1);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k)=semilogy(iter(1:skip:end), err(1:skip:end,k,l,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), err(1:skip:end,k,l,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
% legend(h, {'N_{obs} = 2','N_{obs} = 3', 'N_{obs} = 4'});
xlabel('Iteration n'); ylabel('$|\phi^{(t)}(x) -\phi^{(n)}(x)|_2/|\phi^{(t)}(x)|_2$','interpreter','latex');xlim([0 iter_max]);
ylim([1e-6 1e0]);
 print -depsc2 err_case_3_0_375.eps
hold off;

figure(2);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k)=semilogy(iter(1:skip:end), grad(1:skip:end,k,l,1)/grad(1,k,l,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), grad(1:skip:end,k,l,2)/grad(1,k,l,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$|\nabla J^{(n)}(x)|_2/|\nabla J^{(0)}(x)|_2$','interpreter','latex');xlim([0 iter_max]);
% print -depsc2 err_grad.eps
hold off;

figure(3);
for j_obs = n_obs
    k = find(n_obs==j_obs);
    h(k) = semilogy(iter(1:skip:end), cost(1:skip:end,k,l,1)/cost(1,k,l,1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
           semilogy(iter(1:skip:end), cost(1:skip:end,k,l,2)/cost(1,k,l,2), marker2{k},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$J^{(n)}/J^{(0)}$','interpreter','latex');xlim([0 iter_max]);
print -depsc2 cost_case_3_0_375.eps
hold off;

figure(4);
plot(X,eta_optimum(:,m,l,1),'k-','linewidth',1.5);grid on; hold on;
plot(X,eta_optimum(:,m,l,2),'k--','linewidth',1.5);
set(gca,'FontSize',16);
legend('Linear','Nonlinear');
xlabel('$x$', 'interpreter','latex'); ylabel('$\phi^{(b)}(x)$', 'interpreter', 'latex');xlim([min(X) max(X)]);
print -depsc2 optimum_eta_case_3_fig6e_0_48.eps
hold off;

figure(5); plot(X,eta_exct0(:,m,l,2), 'linewidth',2); hold on
plot(X, (bathymetry(:,m,l,2)-1), 'linewidth', 2);
str1 = '$\eta = \frac{1}{20} exp^{ -(10x)^2}$';
dim1 = [0.32 0.55 0.3 0.3];
str2 = '$\beta = \frac{1}{20}\Big(\tanh\left(2(x+2)\right) - \tanh\left(2(x-2)\right)\Big)$';
dim2 = [0.4 0.000006 0.3 0.3];
annotation('textbox',dim1,'String',str1,'FitBoxToText','on', 'EdgeColor','none', 'interpreter','latex', 'FontSize',24)
annotation('textbox',dim2,'String',str2,'FitBoxToText','on', 'EdgeColor','none', 'interpreter','latex', 'FontSize',24)
xlabel('$x$','interpreter','latex');
hold off;


eta_opt_l = eta_optimum(:,m,l,1);
eta_opt_nl = eta_optimum(:,m,l,2);
eta_exct_spec = eta_exct0(:,m,l,2); 
y_l = abs(fft(eta_opt_l)/length(eta_opt_l));
y_nl = abs(fft(eta_opt_nl)/length(eta_opt_nl));
y_ex = abs(fft(eta_exct_spec)/length(eta_exct_spec));
y_l = y_l(1:length(y_l)/2 +1);
y_nl = y_nl(1:length(y_nl)/2 +1);
y_ex = y_ex(1:length(y_ex)/2 +1);

figure(7);
semilogy(0:50,y_l(1:51),'k-','linewidth',1.5); grid on; hold on;
semilogy(0:50,y_nl(1:51),'k--','linewidth',1.5);
semilogy(0:50,y_ex(1:51),'k-','linewidth',1.5)
hold off;
legend('Linear','Nonlinear','Exact');
xlabel('wavenumber $k$','interpreter', 'latex');
ylabel('Energy Spectrum','interpreter', 'latex');
set(gca,'FontSize',16);
print -depsc2 spec_case_2_fig6f_0_48.eps
