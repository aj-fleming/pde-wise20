N = 1000;
x = linspace(u0(1,1),u0(end,2),N);
figure
plot(x,u(1,:));
xlim([u0(1,1) u0(end,2)])
ylim([0,1])

for idx = 2:size(u,1)
    plot(x,u(idx,:));
    pause(0.001);
    xlim([u0(1,1) u0(end,2)])
    ylim([0,1])
end