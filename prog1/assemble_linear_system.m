% unsure why the program structure is so segmented

%{
assemble_linear_system builds the system Ax=b for the presented problem
it uses linear hat functions to create a basis for V^H_N
therefore the only thing it needs is the grid points to calculate the various
h_i values

assemble_linear_system takes two argument:
  the number of grid points N
  the grid points as a 1xN+1 matrix

assemble_linear_system returns two values:
  the NxN stiffness matrix A
  the Nx1 RHS vector b
%}

function [A, b] = assemble_linear_system(N, grid)
  A = zeros(N, N);
  b = zeros(N, 1);
  for i=1:N
    for j=1:N
      %previously known computation of a(phi_i, phi_i)
      if i == j
        A(j,i) = 1/(grid(i+1)-grid(i)) + 1/(grid(i+2)-grid(i+1));
      
      %previously known computation of a(phi_i, phi_i+1)
      elseif j == i+1
        A(j,i) = -1/(grid(i+2)-grid(i+1));
      
      %previously known computation of a(phi_i, phi_i-1)
      elseif j == i-1
        A(j,i) = -1/(grid(i+1)-grid(i));
      end
    end
    %previously known computation of F(phi_i)
    b(i,1) = 0.5*((grid(i+2)-grid(i+1))+(grid(i+1)-grid(i)));
  end
end