function F_LF = F_LF(uL,uR,lambda)
    F_LF = 0.5*(A(uL)*uL + A(uR)*uR - (uR-uL)/lambda);
end