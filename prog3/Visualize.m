N = 500; % grid points
T = 1; %Time
dt = 0; % time interval; specify it as 0 for the default value to be set based on CFL condition
u0_interval = [-1, 1]; % the domain interval
x_v = linspace(u0_interval(1),u0_interval(2),N); 
% function handle of initial condition pieceiwise function
u0_fun = @(x) 1.*(( x>=-1 && x<-1/2) || ( x>=1/2 && x<=1)); 
% u0_fun = @(x) sin(pi*x); % function handle of initial condition sin(pi x)
F_type = 'LW'; % flux type options: 'naive', 'LF', 'LW'
u = solve_pde(N,T,dt,u0_interval,u0_fun,F_type); % solve the pde
% animated plot
figure
ax1 = axes;
xlim(ax1,[u0_interval(1) u0_interval(2)])
ylim(ax1,[0,1])
size_max = size(u,1);
for idx = 2:size(u,1)
    plot(ax1,x_v,u(idx,:));
    title(ax1,['T = ', num2str(idx/size_max*T)]);
    pause(0.001);
end