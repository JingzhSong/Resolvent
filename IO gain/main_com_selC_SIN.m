%% Calculate the singular values for different Ck and Cd
clear
clc
flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
wall = 'com'; % 'com', 'rigid'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'

kx = 12;kz = 120;c = 10;  % near-wall cycle
%kx = 1;kz = 10;c = 16;   % VLSMs
omega = c*kx;

N = 200; % wall-normal grids
NY = 100; % grids for compliant parameter Cd and Ck
R = 2000;  %Renolds number
i = sqrt(-1);

Cm = 2;
if kx == 12
Ckvector = linspace(27000,31000,NY);
Cdvector = linspace(-10,10,NY);
end
if kx == 1
    Ckvector = linspace(350,700,NY);
    Cdvector = linspace(-10,10,NY);
end


Normcom = zeros(NY,NY);
parpool(50);
parfor Ckk = 1:NY
    Normtmp = zeros(1,NY);
    Ck = Ckvector(Ckk);
    for Cdd = 1:NY
        Cd = Cdvector(Cdd);
        fprintf('Ck=%d,Cd=%d\n',Ckk,Cdd)
        
        ReY = -(omega^2*Cd)/((Ck-omega^2*Cm)^2+(omega*Cd)^2);
        ImY = omega*(Ck-omega^2*Cm)/((Ck-omega^2*Cm)^2+(omega*Cd)^2);  
        Y = ReY + i*ImY;
        
        Normtmp(Cdd) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
    end
    Normcom(Ckk,:) = Normtmp;
end
delete(gcp('nocreate'));
path = [method,'_',eddy,'_',wall,'/singular_CkCd_NY=',num2str(NY),'_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kx),'_kz=',num2str(kz),'_c=',num2str(c)];
save([path,'.mat'],'Normcom');

