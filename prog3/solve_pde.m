function [u , x_vec, t] = solve_pde(N,T,dt,dt_CFL,u0_interval,u0,F_type,f_type)
% Inputs:
%       N - number of gridpoints in space domain
%       T - end time 
%       dt - custom time step; set to 0 if it is not to be defined
%       dt_CFL - CFL factor from 0 to 1;
%       u0_interval - domain of the space e.g: [-1 1]
%       u0 - function handle of the initial condition
%       F_type - numerical flux type; options: 'naive', 'LF', 'LW' 
% Outputs:
%       u - the solution in matrix form; each row is u(x) at specific time
% create the space discretization
x_vec = linspace(u0_interval(1),u0_interval(2),N);
dx = x_vec(2)-x_vec(1);
% set dt according to CFL condition
switch f_type
    case 'advection'
        if dt==0 || dt > dx/2
            dt = dt_CFL*dx/2;
        end
    case 'burgers'
        if dt==0 || dt> dx 
            dt = dt_CFL*dx;
        end
end
% create time discretization vector
t = 0:dt:T;
% initialize u
u = zeros(length(t),N);
% set the initial conditions for u0 (aka u(1,:)) 
for idx2 = 1:length(u(1,:))
    u(1,idx2) = u0(x_vec(idx2));
end
% solve the pde for u at each j and n
for n=1:size(u,1)-1
    for j=1:size(u,2)
        if j > 1 && j < size(u,2)
            % solution for interior space points
            u(n+1,j) = u(n,j) + dt/dx*(F(u(n,j-1),u(n,j),dx,dt,F_type,f_type)-F(u(n,j),u(n,j+1),dx,dt,F_type,f_type));
        elseif j == 1
            % solution for j==1; make it periodic (assign for ul = u(last))
            u(n+1,j) = u(n,j) + dt/dx*(F(u(n,end),u(n,j),dx,dt,F_type,f_type)-F(u(n,j),u(n,j+1),dx,dt,F_type,f_type));
           
        elseif j == size(u,2)
            % solution at the right end of space discretization; (assign for ur = u(1))
            u(n+1,j) = u(n,j) + dt/dx*(F(u(n,j-1),u(n,j),dx,dt,F_type,f_type)-F(u(n,j),u(n,1),dx,dt,F_type,f_type));
        end
    end
end

end
function result = F(ul,ur,dx,dt,F_type,f_type)
% Numerical flux function
% Inputs:
%       ul - u left
%       ur - u right
%       F_type - numerical flux function type; options: 'naive', 'LF', 'LW'
%       f_type - flux function type; options: 'advection', 'burgers'
% Outputs:
%       resuls - scalar value F(ul,ur)
switch F_type
    case 'naive'
        result = 1/2*(f(ul,f_type)+f(ur,f_type));
    case 'LF'
        result = 1/2*(f(ul,f_type)+f(ur,f_type))- dx/2/dt*(ur-ul);
    case 'LW'
        if (ur-ul) == 0 % make sure there is no division by 0
            result = 1/2*(f(ul,f_type)+f(ur,f_type));
        else
            result = 1/2*(f(ul,f_type)+f(ur,f_type))- dt/2/dx*(ur-ul)*((f(ul,f_type)-f(ur,f_type))/(ul-ur))^2;
        end
end
end
function result = f(x,f_type)
% flux function 
switch f_type
    case 'advection'
        result = 2*x;
    case 'burgers'
        result = 1/2*x^2;
end
end