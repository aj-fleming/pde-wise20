function flux = f(u) %u[3x1]
    Cv = 1.4;
%     if u(1)==0 %case divide by 0
%         flux = [u(2);
%             u(3)*(Cv-1)+(3-Cv)/2*(u(2).^2);
%             u(2).*u(3)*Cv-0.5*(Cv-1)*(u(3).^3)];
%     else
        flux = [u(2);
            u(3)*(Cv-1)+(3-Cv)/2*(u(2).^2)./u(1);
            u(2).*u(3)./u(1)*Cv-0.5*(Cv-1)*(u(3).^3)/(u(1).^2)];
%     end
end
