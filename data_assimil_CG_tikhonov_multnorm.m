function [eta_optimum, eta_exct0, beta, X, err, err_std, grad, grad_std, cost, cost_std, n_iter] = ...
    data_assimil_CG_tikhonov_multnorm(nl,N,n_obs,rand_x0,x0_min,dmu,ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, conj_grad_type, iter_chunk, noisy_obs, tikh_par, reg_norm, bath)
% Usage: 
% [eta_optimum, eta_exct0, X, err, err_std, grad, grad_std, cost, cost_std] = ...
%    data_assimil(nl,N,n_obs,rand_x0,x0_min,dmu,ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt)
%
% Performs data assimilation 
%
% Input:
% nl          = nonlinear or linear equations 
% N           =  Number of grid points
% n_obs       = Number of observation points
% rand_x0     =  Use random observation points
% x0_min      = Location of first observation point
% dmu         = Spacing between observation points
% ntrial      = Number of trials
% tmax        =  control time
% xmax        =  domain size: xmin = -xmax
% iter_max    =  Max number of iterations
% line_min    =  Find optimal step
% smooth_grad =  Smooth L2 gradient to Sobolev gradient
% filt        = filtering parameter for Sobolev gradient   
%
% Output
% eta_optimum = final optimized estimate for initial conditions
% eta_exact0  = true initial conditions
% X           = grid point locations
% err         = L2 error for each iteration

% Sample parameters for the data assimilation calculations
% nl          = true;      % nonlinear or linear equations 
% n_obs       = 5;          % Number of observation points
% rand_x0     = false;      % Use random observation points
% x0_min      = 0.2;        % Location of first observation point
% dmu         = 0.2;       % Spacing between observation points
% ntrial      = 1;          % Number of trials
% N           = 512;        % Number of grid points
% tmax        = 2;          % control time
% xmax        = 6;          % domain size: xmin = -xmax
% iter_max        = 500;        % Max number of iterations
% line_min    = true;       % Find optimal step
% smooth_grad = true;      % Smooth L2 gradient to Sobolev gradient
% filt        = 0.4;        % filtering parameter for Sobolev gradient   

% Observation points
x0_max = x0_min+(n_obs-1)*dmu;
x0_vals = x0_min:dmu:x0_max; % Evenly spaced

xmin = -xmax;
h = (xmax - xmin)/N;
dt = h; % Time step (CFL = 1)
tmin = 0;

numSteps = round(abs((tmax-tmin)/dt));
tmax = numSteps*dt;
X = (xmin:h:xmax-h)';

% Set bathymetry
if bath
    b_amp = 0.05;
      beta = b_amp.*exp(-(X-1.5).^2/0.5^2);
%     beta = b_amp.*(tanh(2*(X+2*ones(N,1))) - tanh(2*(X-2*ones(N,1))) );
else
    beta = zeros(N,1);
end


% Set forward and backward solvers for linear and nonlinear cases
if nl
    fw = @fw_SWE_NL_periodic_with_bath;
    bw = @bw_SWE_periodic_mult_NL_v2;
else
    fw = @fw_SWE_periodic;
    bw = @bw_SWE_periodic_mult;
end

for jtrial=1:ntrial
    disp(['Trial ', sprintf('%i', jtrial)]);
    if rand_x0
        dx0_min=0;
        while dx0_min < 2*h
            x0_vals = x0_min + (x0_max-x0_min)*rand(n_obs,1);
            for j=1:n_obs-1
                dx0(j) = x0_vals(j+1)-x0_vals(j);
            end
            dx0_min=min(dx0);
        end
    end
%     x0_vals

    
    % Exact initial conditions
    amp = 0.05; eta_exct0 = amp.*exp(-X.^2/0.1^2);
    [x0_inds, ~] = mult_x0(X,x0_vals);

    
    % STEP I: Run forward solver using exact initial conditions to get observations
    
    u0 = zeros(N,1);
    H = [eta_exct0 ; u0];
    [~, ~, obs_x0_mult, ~, T1, ~] = FW_solve(h, dt, H, 0, tmax, x0_inds, fw, beta);
    % Add gaussian white noise to the observations 
%     mean(obs_x0_mult)
    if noisy_obs
        obs_x0_mult = awgn(obs_x0_mult,20,'measured');
    else
        tikh_par = 0;
    end
    
    % STEP 2: Distort exact IC to get "guess" IC eta0(:,1) used in algorithm with analytical grad_J
    eta0(:,1) = zeros(N,1);
    
    
    % Define eta0(:,2) as initial guess height computed using grad_J
    eta0_2 = zeros(N,1);
    eta0_2(:,1) = eta0(:,1);
    
    % STEP 2.2: Run forward and backwards solvers to get eta_a_t0_x
    H2 = [eta0_2(:,1);u0];
    [eta_all, u_all, eta_x0_t_mult, ~, T2, Y_ex] = FW_solve(h, dt, H2, 0, tmax, x0_inds, fw, beta);
    H3 = zeros(2*N,1);
    [eta_a_t0_x, ~, Y_b] = BW_solve(h, dt, H3, 0, tmax, x0_inds, bw, obs_x0_mult, eta_x0_t_mult, u_all, eta_all, T1, T2, nl);

   
    % Define grad_J
    grad_J = -eta_a_t0_x;
    grad_J = grad_J + tikh_par.*eta0_2(:,1);
    
    % Compute Cost Function
    cost_x = zeros(1,length(T2));
    for i = 1:length(x0_inds)
        cost_x = cost_x + ( (obs_x0_mult(i,:) - eta_x0_t_mult(i,:)).^2) ;
    end
    
    
    if reg_norm == 1    
    reg_term = (norm(eta0_2(:,1))).^2;
    elseif reg_norm == 21
    reg_term = (norm( cent_diff_h(h, eta0_2(:,1) )  )  ).^2;   
    elseif reg_norm == 22
    reg_term = (norm(eta0_2(:,1))).^2 + (norm( cent_diff_h(h, eta0_2(:,1) )  )  ).^2;
    elseif reg_norm == 31
    reg_term = (norm( cent_diff_h(h, cent_diff_h(h, eta0_2(:,1))  )  )  ).^2;
    elseif reg_norm == 32
    reg_term = (norm(eta0_2(:,1))).^2 + (norm( cent_diff_h(h, eta0_2(:,1) )  )  ).^2 + (norm( cent_diff_h(h, cent_diff_h(h, eta0_2(:,1))  )  )  ).^2;
    end
    
    cost(1,jtrial) = trapz(T2,0.5*cost_x) +  (tikh_par/2)*reg_term;
    grad(1,jtrial) = norm(grad_J);
    err(1,jtrial)  = norm(eta_exct0-eta0)/norm(eta_exct0);
    
    % Initial estimate for gradient descent stepsize
    tau_n = 1/length(x0_vals);

    % BEGIN LOOP
    iter=0;
    while norm(grad_J) >= 1e-20 && iter<iter_max
        iter = iter+1;
        % STEP 3
        % Run forward solver to get height at x0 for all t given eta0_1
        H2 = [eta0_2(:,iter) ; u0];
        [eta_all, u_all, eta_x0_t_mult, ~, T2, ~] = FW_solve(h, dt, H2, 0, tmax, x0_inds, fw, beta);
                
        % STEP 4
        % Run backward solver to get adjoint height at t0 for all x
        H3 = zeros(2*N,1);
        [eta_a_t0_x, ~, ~] = BW_solve(h, dt, H3, 0, tmax, x0_inds, bw, obs_x0_mult, eta_x0_t_mult, u_all, eta_all, T1, T2, nl);
        
        % STEP 5: Define gradient of cost function
        grad_J = -eta_a_t0_x;
        grad_J = grad_J + tikh_par.*eta0_2(:,iter);
        steepest_direction(:,iter) = - grad_J;
        
        % STEP 5.5: Smooth gradient
        if smooth_grad
            grad_J = grad_smooth(grad_J,filt);
        end
        
        % for the first iteration, this is steepest descent
        if mod(iter,iter_chunk) == 1 
            % The conjugate direction is the steepest direction at
            % the first iteration
            conj_direction(:,iter) =  steepest_direction(:,iter); 

        
        % STEP 6: Line Minimisation for optimal step size tau_n
        if line_min
            warning('off','optim:fminunc:SwitchingMethod')
            options = optimoptions(@fminunc,'Display', 'off', 'Tolfun',1e-6,'TolX', 1e-6);
            tau0 = tau_n;
            f = @(tau)line_min_mult(obs_x0_mult, eta_a_t0_x, N, h, x0_inds, tau, eta0_2(:,iter), dt, tmin, tmax, fw, beta);
            tau_n = fminunc(f,tau0,options);
            %tau_n = min(tau_n, 3/n_obs);
        end
        
        % STEP 7: Steepest Descent Algorithm
        eta0_2(:,iter+1)= eta0_2(:,iter) + tau_n*conj_direction(:,iter);
        eta_optimum = eta0_2(:,iter+1);
       
        % STEP 8: Run forward solver to get height at x0 for all t given optimised IC
        H4 = [eta0_2(:,iter+1);u0];
        [~, ~, eta_x0_t_mult2, ~, T2, ~] = FW_solve(h, dt, H4, 0, tmax, x0_inds, fw, beta);
                 
        % Compute Cost Function
        cost_x = zeros(1,length(T2));
        for i = 1:length(x0_inds)
            cost_x = cost_x + (obs_x0_mult(i,:) - eta_x0_t_mult2(i,:)).^2 ;
        end
        
        if reg_norm == 1 
        reg_term = (norm(eta0_2(:,iter+1))).^2; 
        elseif reg_norm == 21
        reg_term = (tikh_par/2)*(norm( cent_diff_h(h, eta0_2(:,iter+1))  ) ).^2;   
        elseif reg_norm == 22
        reg_term = (norm(eta0_2(:,iter+1))).^2 + (norm( cent_diff_h(h, eta0_2(:,iter+1) )  )  ).^2;
        elseif reg_norm == 31
        reg_term = (norm( cent_diff_h(h, cent_diff_h(h, eta0_2(:,iter+1))  )  )  ).^2;
        elseif reg_norm == 32
        reg_term = (norm(eta0_2(:,iter+1))).^2 + (norm( cent_diff_h(h, eta0_2(:,iter+1) )  )  ).^2 + (norm( cent_diff_h(h, cent_diff_h(h, eta0_2(:,iter+1))  )  )  ).^2;
        end

        cost(iter+1,jtrial) = trapz(T2,0.5*cost_x)+(tikh_par/2)*reg_term ;
        grad(iter+1,jtrial) = norm(grad_J);
        err(iter+1,jtrial) = norm(eta_exct0-eta_optimum)/norm(eta_exct0);
        disp([sprintf('%1.f',iter),'  ',sprintf('%0.3e',tau_n), ' ', sprintf('%0.3e',norm(grad(iter,jtrial))),'  ', ...
        sprintf('%0.3e',norm(err(iter,jtrial))),' ',sprintf('%0.3e',(tikh_par/2)*reg_term)]);
        
        else
            if conj_grad_type == 0
               % Steepest Descent
               conj_direction(:, iter) = steepest_direction(:,iter);         
        elseif conj_grad_type == 1 
               % Fletcher-Reeves
               b = (steepest_direction(:,iter).' * steepest_direction(:,iter)) / (steepest_direction(:,iter-1).' * steepest_direction(:,iter-1));
               conj_direction(:, iter) = steepest_direction(:,iter) + b*conj_direction(:,iter-1);
        elseif conj_grad_type == 2
%              % Polak-Ribière
               b = (steepest_direction(:,iter).' *  ((steepest_direction(:,iter)) - steepest_direction(:,iter-1) )) / (steepest_direction(:,iter-1).' * steepest_direction(:,iter-1));
               conj_direction(:, iter) = steepest_direction(:,iter) + b*conj_direction(:,iter-1);
        end
                   
        
        
        % STEP 6: Line Minimisation for optimal step size tau_n
        if line_min
            warning('off','optim:fminunc:SwitchingMethod')
            options = optimoptions(@fminunc,'Display', 'off', 'Tolfun',1e-6,'TolX', 1e-6);
            tau0 = tau_n;
            f = @(tau)line_min_mult(obs_x0_mult, eta_a_t0_x, N, h, x0_inds, tau, eta0_2(:,iter), dt, tmin, tmax, fw, beta);
            tau_n = fminunc(f,tau0,options);
            %tau_n = min(tau_n, 3/n_obs);
        end

        % STEP 7: Conjugate gradient method for next optimised guess
        eta0_2(:,iter+1)= eta0_2(:,iter) + tau_n*conj_direction(:,iter);
        
        eta_optimum = eta0_2(:,iter+1);
       
        % STEP 8: Run forward solver to get height at x0 for all t given optimised IC
        H4 = [eta0_2(:,iter+1);u0];
        [~, ~, eta_x0_t_mult2, ~, T2, ~] = FW_solve(h, dt, H4, 0, tmax, x0_inds, fw, beta);
                 
        % Compute Cost Function
        cost_x = zeros(1,length(T2));
        for i = 1:length(x0_inds)
            cost_x = cost_x + (obs_x0_mult(i,:) - eta_x0_t_mult2(i,:)).^2 ;
        end
        
        if reg_norm == 1 
        reg_term = (norm(eta0_2(:,iter+1))).^2;  
        elseif reg_norm == 21
        reg_term = (tikh_par/2)*(norm( cent_diff_h(h, eta0_2(:,iter+1))  ) ).^2;   
        elseif reg_norm == 22
        reg_term = (norm(eta0_2(:,iter+1))).^2 + (norm( cent_diff_h(h, eta0_2(:,iter+1) )  )  ).^2;
        elseif reg_norm == 31
        reg_term = (norm( cent_diff_h(h, cent_diff_h(h, eta0_2(:,iter+1))  )  )  ).^2;
        elseif reg_norm == 32
        reg_term = (norm(eta0_2(:,iter+1))).^2 + (norm( cent_diff_h(h, eta0_2(:,iter+1) )  )  ).^2 + (norm( cent_diff_h(h, cent_diff_h(h, eta0_2(:,iter+1))  )  )  ).^2;
        end
        
        cost(iter+1,jtrial) = trapz(T2,0.5*cost_x) + (tikh_par/2)*reg_term ;
        grad(iter+1,jtrial) = norm(grad_J);
        err(iter+1,jtrial) = norm(eta_exct0-eta_optimum)/norm(eta_exct0);
        disp([sprintf('%1.f',iter),'  ',sprintf('%0.3e',tau_n), ' ', sprintf('%0.3e',norm(grad(iter,jtrial))),'  ', ...
        sprintf('%0.3e',norm(err(iter,jtrial))),' ',sprintf('%0.3e',(tikh_par/2)*reg_term)]);
     end
  end
end

n_iter = iter+1;

err_std  = std(err(end,:));
err      = mean(err,2);

grad_std = std(grad(end,:));
grad     = mean(grad,2);

cost_std = std(cost(end,:));
cost     = mean(cost,2);

disp(['N_obs = ', sprintf('%i', n_obs), '  mean error = ', sprintf('%0.4e',err(end)),' std dev = ', sprintf('%0.4e',err_std(end))]);
