function Au = A(u)
    Cv = 1.4;
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
%     if (u1 == 0) % case divide by 0
%         Au = [0, 1, 0;
%              (Cv-3)/2*(u2)^2, (3-Cv)*(u2/u1), Cv-1;
%              -Cv*u2*u3+(Cv-1)*(u2/u1)^3, Cv*u3/u1-1.5*(Cv-1)*(u2/u1)^2, u2/u1*Cv];
%     else
        Au = [0, 1, 0;
             (Cv-3)/2*(u2/u1)^2, (3-Cv)*(u2), Cv-1;
             -Cv*u2*u3/u1^2+(Cv-1)*(u2)^3, Cv*u3/u1-1.5*(Cv-1)*(u2)^2, u2*Cv];
%     end
end