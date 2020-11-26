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

analytical_1 = @(x) -0.5 .*(x-0.5).^2 + 0.125;
analytical_2 = @(x) (x < 0.5) .* (x) + (x >= 0.5) .* (1-x);

%Set number of hat functions to use
N = 4;
ntries = 100;

width = 0*(1:ntries);
err1s = 0*(1:ntries);
err2s = 0*(1:ntries);
for i=1:ntries
  
  gp = assemble_grid(N, problemtype);
  u_h_1 = zeros(1,length(gp));
  u_h_2 = zeros(1,length(gp));
  
  [A, b] = assemble_linear_system(N, gp);
  b2 = assemble_rhs_dirac(gp, N);

  u_h_1(2:length(gp)-1) = A\b;
  u_h_2(2:length(gp)-1) = A\b2;

  err1s(i) = compute_error(analytical_1, u_h_1, gp);
  err2s(i) = compute_error(analytical_2, u_h_2, gp);
  widths(i) = 1/(N+1);
  N = N + 4;
end

figure(1)
plot(gp, u_h_1);
title("Numerical Solution to (1)");
grid
figure(2)
plot(gp, u_h_2);
title("Numerical Solution to (2)");
figure(3)
loglog(widths, err1s);
title("Error plot for Solution to (1)");
xlabel("grid width");
ylabel("error in the L2 norm")
grid
figure(4)
loglog(widths, err2s);
title("Error plot for Solution to (2)");
xlabel("grid width");
ylabel("error in the L2 norm")
grid
