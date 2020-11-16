% this function finds the RHS of Au=b for problem 2

function [b] = assemble_rhs_dirac(gp, N)
  b = zeros(N, 1);
  for j=1:N
    if 0.5 < gp(j+1) && 0.5 > gp(j)
      b(j,1) = 2*(0.5-gp(j))/(gp(j+1)-gp(j));
    elseif 0.5 < gp(j+2) && 0.5 > gp(j+1)
      b(j,1) = 2*(gp(j+2)-0.5)/(gp(j+2)-gp(j+1));
    end
  end
end
