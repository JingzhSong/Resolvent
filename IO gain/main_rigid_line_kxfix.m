clear
clc

N = 200;  % wall-normal grids
R = 2000; %Renolds number
jz=50; % kz grids

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'
wall = 'rigid';

i = sqrt(-1);

Y=0;

lambdaxVector = [1000,2*R,6*R];
lambdazVector = logspace(0.5,6.5,jz);
singular = zeros(jz,1);
cmax = zeros(jz);

Uc= 24.0126;
Nc = 30;
c = linspace(0,Uc,Nc);

for x = 1:length(lambdaxVector)
    lambdax= lambdaxVector(x);
    lambdax
    for z = 1:jz
        z
        lambdaz = lambdazVector(z);
        kx = 2*pi/lambdax*R;
        
        kz = 2*pi/lambdaz*R;
        singular_vector = zeros(1,Nc);
        for k = 1:Nc
            singular_vector(k) = singularvalue_up(kx,kz,c(k),N,R,flowtype,method,eddy,wall,Y);
        end
        
        [sinTmp,maxTmp] = max(singular_vector);
        cTmp = c(maxTmp);
        singularTmp= sinTmp;
        
        singular(z) = singularTmp;
        cmax(z) =  cTmp;
    end
    
    path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
    save([path,'.mat'],'singular');
    path = [method,'_',eddy,'_',wall,'/cmax_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
    save([path,'.mat'],'cmax');
end



