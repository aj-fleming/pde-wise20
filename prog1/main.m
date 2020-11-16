% Programming Assignment 1
% PDE WiSe 2020
% Submission for: 
% Catalin Dascaliuc 410957
% Alexander J. Fleming 397466
% Sarina Lebs 419937
% Kenny Yung 416070

EQUIDISTANT_GRID = 0;
RANDOM_GRID = 1;

problemtype = 1;

analytical_1 = @(x) -0.5*(x-0.5).^2 + 0.125;

%Set number of hat functions to use
N = 15;

gp = assemble_grid(N, problemtype);


u_h_1 = zeros(1,length(gp));
[A, f] = assemble_linear_system(N, gp);
u_h_1(2:length(gp)-1) = A\f;

figure(1)
plot(gp,u_h_1)
hold on
%fplot(analytical_1, [0 1], '-r')
hold off
grid