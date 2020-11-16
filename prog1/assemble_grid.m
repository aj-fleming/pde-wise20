% unsure why the program structure is so segmented

%{
assemble_grid takes two arguments:
  N: the number of divisions of the region (0,1)
  type:
    0 will generate a linear grid
    1 will generate a random grid 

%}

function [gridpoints] = assemble_grid(N, type)
  if type == 0
    gridpoints = linspace(0, 1, N+2);
  elseif type == 1
    gridpoints = zeros(1, N+2);
    gridpoints(N+2) = 1;
    gridpoints(2:N+1) = sort(rand(1, N));
  end
endfunction
