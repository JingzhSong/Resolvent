clear
clc

N = 200;  % wall-normal grids
R = 2000; %Renolds number
jx=50; % kx grids

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'
wall = 'rigid';

i = sqrt(-1);

Y=0;

lambdazVector = [100,0.6*R];
lambdaxVector = logspace(0.5,6.5,jx);
singular = zeros(jx,1);
cmax = zeros(jx,1);

Uc= 24.0126;
Nc = 30;
c = linspace(0,Uc,Nc);

for z = 1:length(lambdazVector)
    lambdaz = lambdazVector(z);
    lambdaz
    for x = 1:jx
        x
        lambdax = lambdaxVector(x);
        
        kx = 2*pi/lambdax*R;
        kz = 2*pi/lambdaz*R;
        
        singular_vector = zeros(1,Nc);
        for k = 1:Nc
            singular_vector(k) = singularvalue_up(kx,kz,c(k),N,R,flowtype,method,eddy,wall,Y);
        end
        
        [sinTmp,maxTmp] = max(singular_vector);
        cTmp = c(maxTmp);
        singularTmp= sinTmp;
        
        singular(x) = singularTmp;
        cmax(x) =  cTmp;
    end
    
    path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
    save([path,'.mat'],'singular');
    path = [method,'_',eddy,'_',wall,'/cmax_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
    save([path,'.mat'],'cmax');
end



