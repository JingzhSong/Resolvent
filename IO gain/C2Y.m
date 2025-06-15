function [ReY,ImY] = C2Y(Cm,Ck,Cd,omega)

ReY = -(omega.^2.*Cd)./((Ck-omega.^2.*Cm).^2+(omega.*Cd).^2);
ImY = omega.*(Ck-omega.^2.*Cm)./((Ck-omega.^2.*Cm).^2+(omega.*Cd).^2);

% fprintf('ReY=%d,ImY=%d',ReY,ImY);
end

