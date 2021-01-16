function u = solve_pde(N,T,u0,F_type)
% u0 - matrix size: n by 3; n number of intervals; 1st value left bound,
% 2nd value right bound, 3rd value- value of function2
f_type = 'burgers';
% initialize
interval = linspace(u0(1,1),u0(end,2),N);
dx = interval(2)-interval(1);
switch f_type
    case 'advection'
        dt = dx/2;
    case 'burgers'
        dt = dx/u0(end,2);
end
t = 0:dt:T;
u = zeros(length(t),N);
for idx = 1:size(u0,1)
    for idx2 = 1:length(u(1,:))
        if interval(idx2)>=u0(idx,1) && interval(idx2)<u0(idx,2)
            u(1,idx2) = u0(idx,3);
        end
    end
end
u(1,end) = 1;
for n=1:size(u,1)
    for j=1:size(u,2)
        if j > 1 && j < size(u,2)
            u(n+1,j) = u(n,j) + dt/dx*(F(u(n,j-1),u(n,j),dx,dt,F_type,f_type)-F(u(n,j),u(n,j+1),dx,dt,F_type,f_type));
        elseif j == 1
            u(n+1,j) = u(n,j) + dt/dx*(F(u(n,end),u(n,j),dx,dt,F_type,f_type)-F(u(n,j),u(n,j+1),dx,dt,F_type,f_type));
           
        elseif j == size(u,2)
            u(n+1,j) = u(n,j) + dt/dx*(F(u(n,j-1),u(n,j),dx,dt,F_type,f_type)-F(u(n,j),u(n,1),dx,dt,F_type,f_type));
        end
    end
end

end
function result = F(ul,ur,dx,dt,F_type,f_type)
switch F_type
    case 'naive'
        result = 1/2*(f(ul,f_type)+f(ur,f_type));
    case 'LF'
        result = 1/2*(f(ul,f_type)+f(ur,f_type))- dx/2/dt*(ur-ul);
    case 'LW'
        if (ur-ul) == 0
            result = 1/2*(f(ul,f_type)+f(ur,f_type));
        else
            result = 1/2*(f(ul,f_type)+f(ur,f_type))- dt/2/dx*(ur-ul)*((f(ul,f_type)-f(ur,f_type))/(ul-ur))^2;
        end
end
end
function result = f(x,f_type)
switch f_type
    case 'advection'
        result = 2*x;
    case 'burgers'
        result = 1/2*x^2;
end
end