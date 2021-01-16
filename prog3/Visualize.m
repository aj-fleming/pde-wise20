N = 1000;
T = 5;
u0_interval = [-1, 1];
x_v = linspace(u0_interval(1),u0_interval(2),N);
% fun = @(x) 1.*(( x>=-1 && x<-1/2) || ( x>=1/2 && x<1));
fun = @(x) sin(pi*x);
u = solve_pde(N,T,u0_interval,fun,'LF');
figure
plot(x_v,u(1,:));
xlim([u0_interval(1) u0_interval(2)])
ylim([0,1])

for idx = 2:size(u,1)
    plot(x_v,u(idx,:));
    pause(0.001);
    xlim([u0_interval(1) u0_interval(2)])
    ylim([0,1])
end