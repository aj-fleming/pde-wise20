function [U] = Euler_Solver(u0)
    U = zeros(3,N,length(T));
    U(1,:,1) = u0(1,:);
    U(2,i,1) = u0(2,:);
    U(3,i,1) = u0(3,:);
    
end