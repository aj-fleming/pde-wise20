clear all
T = 0.5; %Time
u0 = @(x) 1*(( x>=-3-1/2 && x<-2-1/2) || ( x>=-1-1/2 && x<-1/2) || ( x>=1/2 && x<=1+1/2) || ( x>=2+1/2 && x<=3+1/2));
k = 1;
f_type = 'advection'; % options: 'advection', 'burgers'
dt_CFL = 1; % dt_CFL - CFL factor from 0 to 1;
u0_interval = [-1, 1]; % the domain interval
uEx = @(x,t) u0(x-2*t);
% function handle of initial condition pieceiwise function
for pow = 6:12
    dt = 0;
    N = 2^pow;    
    [u_dt_naive, ~, ~] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,'naive',f_type); % solve the pde
    [u_dt_LF, ~, ~] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,'LF',f_type); % solve the pde
    [u_dt_LW, x_v, t] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,'LW',f_type); % solve the pde
    u = zeros(length(t),length(x_v));
    for idx = 1:length(t)
        for idx2 = 1:length(x_v)
            u(idx,idx2) = uEx(x_v(idx2),(t(idx)));
        end
        u_err_naive(idx) = norm(u_dt_naive(idx,:) - u(idx,:),1);
        u_err_LF(idx) = norm(u_dt_LF(idx,:) - u(idx,:),1);
        u_err_LW(idx) = norm(u_dt_LW(idx,:) - u(idx,:),1);
    end
    u_err_pow_naive(k) = max(u_err_naive);
    u_err_pow_LF(k) = max(u_err_LF);
    u_err_pow_LW(k) = max(u_err_LW);
    k = k+1;
end
loglog(2.^(6:12),u_err_pow_naive)
hold on
loglog(2.^(6:12),u_err_pow_LF)
loglog(2.^(6:12),u_err_pow_LW)
legend('naive', 'LF', 'LW')
xlabel('Number of grid points')
ylabel('Error L_1 norm')
save('conv_disc_adv')