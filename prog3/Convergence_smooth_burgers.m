clear all
T = 0.06; %Time
u0 = @(x) sin(pi*x); % function handle of initial condition sin(pi x)
u0_neg = @(x) -sin(pi*x);
a = -1;
b = 1;
[xmin,u0min] = fminbnd(@(x) u0(x),a,b);
[xmax,u0max] = fminbnd(@(x) -u0(x),a,b);
u0max = -u0max;
uEx =@(x,t) fminbnd(@(u) (u - u0(x-u*t)).^2,u0min,u0max);
k = 1;
f_type = 'burgers'; % options: 'advection', 'burgers'
dt_CFL = 1; % dt_CFL - CFL factor from 0 to 1;
u0_interval = [-1, 1]; % the domain interval
% function handle of initial condition pieceiwise function
% u0_fun = @(x) 1*(( x>=-1 && x<-1/2) || ( x>=1/2 && x<=1));
u0_fun = @(x) sin(pi*x); % function handle of initial condition sin(pi x)
for pow = 6:12
    dt = 0;
    N = 2^pow;    
    [u_dt_naive, ~, ~] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,'naive',f_type); % solve the pde
    [u_dt_LF, ~, ~] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,'LF',f_type); % solve the pde
    [u_dt_LW, x_v, t] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,'LW',f_type); % solve the pde
    u = zeros(length(t),length(x_v));
    dx = x_v(2)-x_v(1);
    dt = t(2) - t(1);
    u_p0 = zeros(length(t),length(x_v));
    max_xv = length(x_v);
    for idx = 1:length(t)
        for idx2 = 1:max_xv
            u(idx,idx2) = uEx(x_v(idx2),t(idx));
            if idx2 ==1
                u_p0(idx,idx2) = uEx(x_v(idx2),t(idx));
            elseif idx2==max_xv
                u_p0(idx,idx2) = uEx(x_v(idx2),t(idx));
            else
                u_p0(idx,idx2) = (uEx(x_v(idx2)-dx/2,t(idx)) + uEx(x_v(idx2)+dx/2,t(idx)))/2;
            end
        end
    end
    u_err_pow_naive(k) = norm(u_dt_naive-u,Inf);
    u_err_pow_LF(k) = norm(u_dt_LF-u,Inf);
    u_err_pow_LW(k) = norm(u_dt_LW-u,Inf);
    u_err_pow_naive_p0(k) = norm(u_dt_naive-u_p0,Inf);
    u_err_pow_LF_p0(k) = norm(u_dt_LF-u_p0,Inf);
    u_err_pow_LW_p0(k) = norm(u_dt_LW-u_p0,Inf);
    k = k+1;
    dx_vec(k) = dx;
end
figure
loglog(2.^(6:12),u_err_pow_naive)
hold on
loglog(2.^(6:12),u_err_pow_LF)
loglog(2.^(6:12),u_err_pow_LW)
legend('naive', 'LF', 'LW')
xlabel('Number of grid points')
ylabel('Error L_1 norm')
title('Exact Solution Error')
figure
loglog(2.^(6:12),u_err_pow_naive_p0)
hold on
loglog(2.^(6:12),u_err_pow_LF_p0)
loglog(2.^(6:12),u_err_pow_LW_p0)
legend('naive', 'LF', 'LW')
xlabel('Number of grid points')
ylabel('Error L_1 norm')
title('P0 Error')
save('conv_smooth_burg')