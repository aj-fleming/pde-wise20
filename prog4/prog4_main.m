N = 100; %grid cells
x_start = 0;
x_end = 1;
T_end = 0.5;
dx = (x_end-x_start)/N;
cfl_a = 0.75;
Cv = 1.4;
x = linspace(x_start,x_end,N);
u0 = zeros(3,N); %rho,v,P
IC = '2'; %initial condition
switch IC
    case '1'
        u0 = [1, 0, 1; %rho, v, p
              0.125, 0, 0.1];
    case '2'
        u0 = [1, -1, 0.4;
              1, 1, 0.4];
end
%functions for IC in U[rho; rho*v, E]
U0 = cell(3,1);
U0{1} = @(x) u0(1,1).* and(x>=0, x<=0.5) + u0(2,1).* and(x<=1, x>=0.5);
U0{2} = @(x) u0(1,1)* u0(1,2).* and(x>=0, x<=0.5) + u0(2,1)*u0(2,2).* and(x<=1, x>=0.5);
U0{3} = @(x) 0.5.*U0{1}(x).*U0{2}(x).^2 + u0(1,3)/(Cv-1).* and(x>=0, x<=0.5) + u0(2,3)/(Cv-1).* and(x<=1, x>=0.5);
%calculate cell mean values for U[t=0]
U = zeros(3,N,2);
for i = 1:N
    U(1,i,1) =  integral(U0{1},x_start+(i-1)*dx, x_start+(i*dx))/dx;
    U(2,i,1) =  integral(U0{2},x_start+(i-1)*dx, x_start+(i*dx))/dx;
    U(3,i,1) =  integral(U0{3},x_start+(i-1)*dx, x_start+(i*dx))/dx;
end

time = 0; n = 1; T = [0];
while(true) %iterate over time
    eigs = zeros(1,N);
    for i = 1:N
        eigs(i) = max(abs(fast_eigs(U(:,i,n), Cv)));
    end
    dt = cfl_a*dx/max(eigs); %CFL condition
    if time > T_end
        break;
    end
    time = time + dt;
    T = [T time];
    lambda = dt/dx;
    %calc first one
    U(:,1,n+1) = U(:,1,n)-lambda*(F_HLL(U(:,1,n),U(:,2,n),Cv)-F_HLL(U(:,1,n),U(:,1,n),Cv));
    %calc 2 thru N-1
    for j = 2:size(U,2)-1 %iterate over x
        F_1 = F_HLL(U(:,j,n),U(:,j+1,n),Cv);
        F_2 = F_HLL(U(:,j-1,n),U(:,j,n),Cv);
        U(:,j,n+1) = U(:,j,n)-lambda*(F_1-F_2);
    end
    %calc last one
    j = size(U,2);
    U(:,j,n+1) = U(:,j,n)-lambda*(F_HLL(U(:,j,n),U(:,j,n),Cv)-F_HLL(U(:,j-1,n),U(:,j,n),Cv));
    n = n+1;
end

figure(1)
set(gcf,'position',[100 100 1000 600]);
for t = 1:length(T)
    subplot(3,1,1);
    plot(x,U(1,:,t));
    ylim([0,1]);
    xlabel("T = "+T(t));
    title("Density");
    subplot(3,1,2)
    plot(x,U(2,:,t));
    xlabel("T = "+T(t));
    title("Velocity");
    subplot(3,1,3)
    plot(x,U(3,:,t));
    %ylim([0,1.5]);
    xlabel("T = "+T(t));
    title("Pressure");
    pause(0.0005);
end

function e = fast_eigs(u, heat_ratio)
  v = u(2)/u(1);
  P = (heat_ratio-1)*(u(3)-0.5*u(1)*v^2);
  mach = sqrt(heat_ratio*u(1)/P);
  e = [v-mach, v, v+mach];
end

% lax-friedrichs "cheat" flux to test/compare
% takes two 3x1 arrays and the courant number
%function F_LF = F_LF(uL,uR,lambda)
%    F_LF = 0.5*(A(uL)*uL + A(uR)*uR - (uR-uL)/lambda);
%end

% HLL flux as defined in the lecture
% accepts two 3x1 vectors and the heat capacity ratio (>1)
function FLF = F_HLL(uL,uR,heat_ratio) %u[3x1]
    FLF = zeros(3,1);
    evals = [fast_eigs(uL, heat_ratio) fast_eigs(uR, heat_ratio)];
    aL = min(evals);
    aR = max(abs(evals));
    if aL >= 0
        FLF = f(uL,heat_ratio);
    elseif (aL <= 0) & (aR >= 0)
        FLF = 1./(aR-aL).*(aR.*f(uL,heat_ratio)-aL.*f(uR,heat_ratio)+aL.*aR.*(uR-uL));
    elseif aR >= 0
        FLF = f(uR,heat_ratio);
    end
end
%calc jacobian
function Au = A(u, heat_ratio)
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
    Au = [0, 1, 0;
          (heat_ratio-3)/2*(u2/u1)^2, (3-heat_ratio)*(u2/u1), heat_ratio-1;
          -heat_ratio*u2*u3/(u1^2)+(heat_ratio-1)*(u2/u1)^3, heat_ratio*u3/u1-1.5*(heat_ratio-1)*(u2/u1)^2, u2/u1*heat_ratio];
end
%calc system F(U)
function flux = f(u, heat_ratio) %u[3x1]
  flux = [u(2);
          u(3)*(heat_ratio-1)+(3-heat_ratio)/2*(u(2).^2)./u(1);
          u(2)*u(3)/u(1)*heat_ratio-0.5*(heat_ratio-1)*(u(3).^3)/(u(1).^2)];
end
