function [vel_u,vel_v,vel_w,p,f1,f2,f3,RSy,RS] = UF_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y)
% (kx,kz,c): wave number
% N: y+ grid number
% R: Reynolds number
% flowtype: 'turbu', 'couette', 'poiseu'
% method: 'IOA' or 'SIOA'
% wall: 'rigid' or 'com'

i = sqrt(-1);

[y,DM] = chebdif(N,2);
D1 = DM(1:N,1:N,1);
D2 = DM(1:N,1:N,2);
I = eye(N);
grd = [i*kx.*I; D1; i*kz.*I];
lplc = D2 - kx.^2.*I - kz.^2.*I;

%% velocity profile
if strcmp(flowtype, 'couette')
    U = @(y)y;
    U0 = U(y);
    U1 = D1*U(y);
end

if strcmp(flowtype, 'poiseu')
    U = @(y)1-y.^2;
    U0 = U(y);
    U1 = D1*U(y);
end

if strcmp(flowtype, 'turbu')
    k=0.426;
    alpha=25.4;
    NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
    nuT0 = NuT(y);
    nuT1 = D1*nuT0;
    DUDy = @(y) R.*(-y)./NuT(y);
    U1 = DUDy(y);
    U0 = zeros(N,1);

    for j=1:N
        U0(j) = integral(DUDy,-1,y(j));
    end
end


%% resolvent matrix

if strcmp(eddy,'eddyoff')
    A11 = -1i*kx.*diag(U0) + lplc./R;
    A = [A11,      -diag(U1),  zeros(N),  -1i*kx.*I;
        zeros(N),  A11,        zeros(N),  -D1;
        zeros(N),  zeros(N),   A11,       -1i*kz.*I;
        1i*kx.*I,   D1,        1i*kz.*I,  zeros(N)];
end

if strcmp(eddy,'eddyon')
    A11 = -i*kx*diag(U0) + diag(nuT0)*lplc./R + diag(nuT1)*D1./R;

    A = [A11,     -diag(U1)+i.*kx.*diag(nuT1)./R, zeros(N), -i*kx.*I;
        zeros(N), A11+diag(nuT1)*D1./R,           zeros(N), -D1;
        zeros(N), i.*kz.*diag(nuT1)./R,           A11,      -i*kz.*I;
        i*kx.*I,  D1,                             i*kz.*I,  zeros(N)];
end

if strcmp(method,'IOA')
    C = [eye(3*N),zeros(3*N,N)];
    [~,w]=clencurt(N+1);
    IWC = blkdiag(diag(w(2:N+1).^0.5),diag(w(2:N+1).^0.5),diag(w(2:N+1).^0.5));
    C = IWC*C;
end

if strcmp(method,'SIOA')
    C = [blkdiag(grd,grd,grd),zeros(9*N,N)];
    [~,w]=clencurt(N+1);
    IWC = blkdiag(diag(w(2:N+1).^0.5),diag(w(2:N+1).^0.5),diag(w(2:N+1).^0.5));
    C = blkdiag(IWC,IWC,IWC)*C;
end

B0 = [eye(3*N); zeros(N,3*N)];
IWB = blkdiag(diag(w(2:N+1).^(-0.5)),diag(w(2:N+1).^(-0.5)),diag(w(2:N+1).^(-0.5)));
B = B0*IWB;
E = blkdiag(eye(3*N), zeros(N));


%% boundary conditions

if strcmp(wall,'rigid')
    % u(0)=u(N)=0
    A(1,:)=[1,zeros(1,4*N-1)];
    A(N,:)=[zeros(1,N-1),1,zeros(1,3*N)];
    B(1,:)=zeros(1,3*N);
    B(N,:)=zeros(1,3*N);
    E(1,:)=zeros(1,4*N);
    E(N,:)=zeros(1,4*N);
    % v(0)=v(N)=0
    A(N+1,:)=[zeros(1,N),1,zeros(1,3*N-1)];
    A(2*N,:)=[zeros(1,2*N-1),1,zeros(1,2*N)];
    B(N+1,:)=zeros(1,3*N);
    B(2*N,:)=zeros(1,3*N);
    E(N+1,:)=zeros(1,4*N);
    E(2*N,:)=zeros(1,4*N);
    % w(0)=w(N)=0
    A(2*N+1,:)=[zeros(1,2*N),1,zeros(1,2*N-1)];
    A(3*N,:)=[zeros(1,3*N-1),1,zeros(1,N)];
    B(2*N+1,:)=zeros(1,3*N);
    B(3*N,:)=zeros(1,3*N);
    E(2*N+1,:)=zeros(1,4*N);
    E(3*N,:)=zeros(1,4*N);
end

if strcmp(wall,'com')
    % compliant boundary conditions
    % icu(N)+UyN*Y*p(N)=0;icu(0)-Uy0*Y*p(0)=0;
    A(1,:) = [zeros(1,3*N),Y*DUDy(1),zeros(1,N-1)];
    A(N,:) = [zeros(1,4*N-1),-Y*DUDy(-1)];
    B(1,:) = zeros(1,3*N);
    B(N,:) = zeros(1,3*N);
    
    % v(N)+Yp(N)=0;v(0)-Yp(0)=0;
    A(N+1,:) = [zeros(1,N),1,zeros(1,2*N-1),Y,zeros(1,N-1)];
    A(2*N,:) = [zeros(1,2*N-1),1,zeros(1,2*N-1),-Y];
    B(N+1,:) = zeros(1,3*N);
    B(2*N,:) = zeros(1,3*N);
    E(N+1,:) = zeros(1,4*N);
    E(2*N,:) = zeros(1,4*N);

    % w(0)=w(N)=0
    A(2*N+1,:) = [zeros(1,2*N),1,zeros(1,2*N-1)];
    A(3*N,:) = [zeros(1,3*N-1),1,zeros(1,N)];
    B(2*N+1,:) = zeros(1,3*N);
    B(3*N,:) = zeros(1,3*N);
    E(2*N+1,:) = zeros(1,4*N);
    E(3*N,:) = zeros(1,4*N);
end

%% compute U&F
if strcmp(method,'IOA')
    
    H = C*((-i*c*kx*E-A)\B);
    [U,D,V] = svd(H);
    eig = D(1,1);
    Ho = (-i*c*kx*E-A)\B0;
    
    vel_u = inv(diag(w(2:N+1).^0.5))*U(1:N,1);
    vel_v = inv(diag(w(2:N+1).^0.5))*U(N+1:2*N,1);
    vel_w = inv(diag(w(2:N+1).^0.5))*U(2*N+1:3*N,1);
    
    f1 = inv(diag(w(2:N+1).^0.5))* V(1:N,1);
    f2 = inv(diag(w(2:N+1).^0.5))* V(N+1:2*N,1);
    f3 = inv(diag(w(2:N+1).^0.5))* V(2*N+1:3*N,1);

    f = [f1;f2;f3];
    
    Ho = (-i*c*kx*E-A)\B0;
    uvwp = Ho * f;
    p = uvwp(3*N+1:4*N,1) ./ eig;

    RSy = real(conj(vel_u).*vel_v);
    RSy_eig = RSy.*eig.^2;
    RS = sum(RSy.*w(2:N+1)'.*y).*eig.^2;
end

if strcmp(method,'SIOA')
    
    BlockStructure = [N,3*N;N,3*N;N,3*N];
    omega = c*kx;
    Mp = C*((1i*omega*E-A)\B);
    [bound,muinfo] = mussv(Mp,BlockStructure,'Uf');
    eig = bound(1,1);
    [~,VSigma,~] = mussvextract(muinfo);
    Dl = VSigma.DLeft;
    Dr = VSigma.DRight;
    H = Dl*(C*((1i*omega*E-A)\B))*inv(Dr); 
    [U,~,V] = svd(H);
    vel_u = inv(diag(w(2:N+1).^0.5))*U(1:N,1)./i./kx;
    vel_v = inv(diag(w(2:N+1).^0.5))*U(3*N+1:4*N,1)./i./kx;
    vel_w = inv(diag(w(2:N+1).^0.5))*U(6*N+1:7*N,1)./i./kx;
  
    f1 = inv(diag(w(2:N+1).^0.5))* V(1:N,1);
    f2 = inv(diag(w(2:N+1).^0.5))* V(N+1:2*N,1);
    f3 = inv(diag(w(2:N+1).^0.5))* V(2*N+1:3*N,1);

    f = [f1;f2;f3];
    
    Hoo = ((1i*omega*E-A)\B)*inv(Dr);
    uvwp = Hoo * V(:,1);
    p = uvwp(3*N+1:4*N,1) ./ eig .*Dl(3*N+1,3*N+1); % The same scale as v

    RSy = real(conj(vel_u).*vel_v);
%     RSy_eig = RSy.*eig.^2;
    RS = sum(RSy.*w(2:N+1)'.*y).*eig.^2;
   
end

end

