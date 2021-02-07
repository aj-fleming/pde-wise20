function FLF = F(uL,uR) %u[3x1]
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