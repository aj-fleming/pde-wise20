N = 500; % grid points
T = 0.3; %Time
dt = 0;
dt_CFL = 1; % dt_CFL - CFL factor from 0 to 1;
u0_interval = [-1, 1]; % the domain interval
% function handle of initial condition pieceiwise function
% u0_fun = @(x) 1*(( x>=-1 && x<-1/2) || ( x>=1/2 && x<=1)); 
u0_fun = @(x) sin(pi*x); % function handle of initial condition sin(pi x)
F_type = 'LF'; % flux type options: 'naive', 'LF', 'LW'
f_type = 'burgers'; % options: 'advection', 'burgers'
[u, x_v, t] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0_fun,F_type,f_type); % solve the pde
% animated plot
figure
ax1 = axes;
xlim(ax1,[u0_interval(1) u0_interval(2)])
ylim(ax1,[0,1])
size_max = size(u,1);
for idx = 2:size(u,1)
    plot(ax1,x_v,u(idx,:));
    title(ax1,['T = ', num2str(t(idx))]);
    pause(0.005);
end