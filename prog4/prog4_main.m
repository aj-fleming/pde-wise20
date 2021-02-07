N = 100;
x_start = 0;
x_end = 1;
T_end = 0.25;
dx = (x_end-x_start)/N;
cfl_a = 0.5;
Cv = 1.4;
dt = cfl_a*dx; %/max(abs(eig(Aj))); %CFL condition
lambda = dt/dx;
T = 0:dt:T_end;
if T_end-T(end) > 0
    T = [T, T_end];
end
x = linspace(x_start,x_end,N);
u0 = zeros(3,N); %rho,v,P
IC = '1'; %initial condition
switch IC
    case '1'
        u0 = [1, 0, 1; %rho, v, p
              0.125, 0, 0.1];
        
    case '2'
        u0 = [1, -1, 0.4;
              1, 1, 0.4];
end
U0 = cell(3,1);
U0{1} = @(x) u0(1,1).* and(x>=0, x<=0.5) + u0(2,1).* and(x<=1, x>=0.5);
U0{2} = @(x) u0(1,1)* u0(1,2).* and(x>=0, x<=0.5) + u0(2,1)*u0(2,2).* and(x<=1, x>=0.5);
U0{3} = @(x) 0.5.*U0{1}(x).*U0{2}(x).^2 + u0(1,3)/(Cv-1).* and(x>=0, x<=0.5) + u0(2,3)/(Cv-1).* and(x<=1, x>=0.5);
U = zeros(3,N,length(T));
%calc cell mean value
for i = 1:N
    U(1,i,1) =  quadcc(U0{1},x_start+(i-1)*dx, x_start+(i*dx))/dx;
    U(2,i,1) =  quadcc(U0{2},x_start+(i-1)*dx, x_start+(i*dx))/dx;
    U(3,i,1) =  quadcc(U0{3},x_start+(i-1)*dx, x_start+(i*dx))/dx;
end

for n = 1:size(U,3)-1 %iterate over time
    %calc first one
    U(:,1,n+1) = U(:,1,n)-dt/dx*(F(U(:,1,n),U(:,2,n))-F(U(:,1,n),U(:,1,n)));
    for j = 2:size(U,2)-1 %iterate over x
        F_1 = F(U(:,j,n),U(:,j+1,n));
        F_2 = F(U(:,j-1,n),U(:,j,n));
        U(:,j,n+1) = U(:,j,n)-dt/dx*(F_1-F_2);
    end
    %calc last one
    j = size(U,2);
    U(:,j,n+1) = U(:,j,n)-dt/dx*(F(U(:,j,n),U(:,j,n))-F(U(:,j-1,n),U(:,j,n)));
end





