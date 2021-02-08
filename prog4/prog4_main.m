N = 50; %grid cells
x_start = 0;
x_end = 1;
T_end = 0.5;
dx = (x_end-x_start)/N;
cfl_a = 0.5;
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
U0 = cell(3,1);
U0{1} = @(x) u0(1,1).* and(x>=0, x<=0.5) + u0(2,1).* and(x<=1, x>=0.5);
U0{2} = @(x) u0(1,1)* u0(1,2).* and(x>=0, x<=0.5) + u0(2,1)*u0(2,2).* and(x<=1, x>=0.5);
U0{3} = @(x) 0.5.*U0{1}(x).*U0{2}(x).^2 + u0(1,3)/(Cv-1).* and(x>=0, x<=0.5) + u0(2,3)/(Cv-1).* and(x<=1, x>=0.5);
U = zeros(3,N,2);
%calc cell mean value
for i = 1:N
    U(1,i,1) =  integral(U0{1},x_start+(i-1)*dx, x_start+(i*dx))/dx;
    U(2,i,1) =  integral(U0{2},x_start+(i-1)*dx, x_start+(i*dx))/dx;
    U(3,i,1) =  integral(U0{3},x_start+(i-1)*dx, x_start+(i*dx))/dx;
end

time = 0; n = 1; T = [0];
while(true) %iterate over time
    eigs = zeros(1,N);
    for i = 1:N
        Ai = A(U(:,i,n));
        eigs(i) = max(abs(eig(Ai)));
    end
    dt = cfl_a*dx/max(eigs); %CFL condition
    if time > T_end
        break;
    end
    time = time + dt;
    T = [T time];
    lambda = dt/dx;
    %calc first one
    U(:,1,n+1) = U(:,1,n)-lambda*(F(U(:,1,n),U(:,2,n),lambda)-F(U(:,1,n),U(:,1,n),lambda));
    %calc 2 thru N-1
    for j = 2:size(U,2)-1 %iterate over x
        F_1 = F(U(:,j,n),U(:,j+1,n),lambda);
        F_2 = F(U(:,j-1,n),U(:,j,n),lambda);
        U(:,j,n+1) = U(:,j,n)-lambda*(F_1-F_2);
    end
    %calc last one
    j = size(U,2);
    U(:,j,n+1) = U(:,j,n)-lambda*(F(U(:,j,n),U(:,j,n),lambda)-F(U(:,j-1,n),U(:,j,n),lambda));
    n = n+1;
end
dt = T_end-time;
if dt > 0
    %calc last time step
    T = [T T_end];
    lambda = dt/dx;
    U(:,1,n+1) = U(:,1,n)-dt/dx*(F(U(:,1,n),U(:,2,n),lambda)-F(U(:,1,n),U(:,1,n),lambda));
    %calc 2 thru N-1
    for j = 2:size(U,2)-1 %iterate over x
        F_1 = F(U(:,j,n),U(:,j+1,n),lambda);
        F_2 = F(U(:,j-1,n),U(:,j,n),lambda);
        U(:,j,n+1) = U(:,j,n)-dt/dx*(F_1-F_2);
    end
    %calc last one
    j = size(U,2);
    U(:,j,n+1) = U(:,j,n)-dt/dx*(F(U(:,j,n),U(:,j,n),lambda)-F(U(:,j-1,n),U(:,j,n),lambda));
end

figure(1)
set(gcf,'position',[100 100 1000 600]);
for t = 1:length(T)
    tiledlayout(3,1);
    nexttile
    plot(x,U(1,:,t));
    xlabel("T = "+t*dt);
    title("Density");
    nexttile
    plot(x,U(2,:,t));
    xlabel("T = "+t*dt);
    title("Velocity");
    nexttile
    plot(x,U(3,:,t));
    xlabel("T = "+t*dt);
    title("Pressure");
    pause(0.005);
end

function F_LF = F_LF(uL,uR,lambda)
    F_LF = 0.5*(A(uL)*uL + A(uR)*uR - (uR-uL)/lambda);
end
function FLF = F(uL,uR,lambda) %u[3x1]
    FLF = zeros(3,1);
    aL = min(eig(A(uL)),eig(A(uR)));
    aR = max(abs(eig(A(uL))),abs(eig(A(uR))));
    if aL >= 0
        FLF = f(uL);
    elseif (aL <= 0) & (aR >= 0)
        FLF = 1./(aR-aL).*(aR.*f(uL)-aL.*f(uR)+aL.*aR.*(uR-uL));
    elseif aR >= 0
        FLF = f(uR);
    end
end

function Au = A(u)
    Cv = 1.4;
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
%     if (u1 == 0) % case divide by 0
%         Au = [0, 1, 0;
%              (Cv-3)/2*(u2)^2, (3-Cv)*(u2), Cv-1;
%              -Cv*u2*u3+(Cv-1)*(u2)^3, Cv*u3-1.5*(Cv-1)*(u2)^2, u2*Cv];
%     else
        Au = [0, 1, 0;
             (Cv-3)/2*(u2/u1)^2, (3-Cv)*(u2/u1), Cv-1;
             -Cv*u2*u3/(u1^2)+(Cv-1)*(u2/u1)^3, Cv*u3/u1-1.5*(Cv-1)*(u2/u1)^2, u2/u1*Cv];
%     end
end

function flux = f(u) %u[3x1]
    Cv = 1.4;
%     if u(1)==0 %case divide by 0
%         flux = [u(2);
%             u(3)*(Cv-1)+(3-Cv)/2*(u(2).^2);
%             u(2).*u(3)*Cv-0.5*(Cv-1)*(u(3).^3)];
%     else
        flux = [u(2);
            u(3)*(Cv-1)+(3-Cv)/2*(u(2).^2)./u(1);
            u(2)*u(3)/u(1)*Cv-0.5*(Cv-1)*(u(3).^3)/(u(1).^2)];
%     end
end
