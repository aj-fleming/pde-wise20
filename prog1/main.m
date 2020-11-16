% Programming Assignment 1
% PDE WiSe 2020
% Submission for: 
% Catalin Dascaliuc 410957
% Alexander J. Fleming 397466
% Sarina Lebs 419937
% Kenny Yung 416070

EQUIDISTANT_GRID = 0;
RANDOM_GRID = 1;

problemtype = 0;

analytical_1 = @(x) -0.5*(x-0.5).^2 + 0.125;

%Set number of hat functions to use
N = 100;

gp = assemble_grid(N, problemtype);


u_h_1 = zeros(1,length(gp));
u_h_2 = zeros(1,length(gp));

[A, b] = assemble_linear_system(N, gp);
b2 = assemble_rhs_dirac(gp, N);

u_h_1(2:length(gp)-1) = A\b;
u_h_2(2:length(gp)-1) = A\b2;

figure(1)
fplot(analytical_1, [0 1], '-r')
hold on
plot(gp,u_h_1)
hold off
grid

figure(2)
plot(gp,u_h_2)
grid