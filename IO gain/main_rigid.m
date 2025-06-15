clear
clc

N = 200;  % wall-normal grids
R = 2000;  %Renolds number
ix=50; 
jz=50; 
nt=ix*jz;

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
wall = 'com'; % 'com', 'rigid'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'


i = sqrt(-1);

Y=0;

singular = zeros(ix,jz);
cmax = zeros(ix,jz);

Uc= 24.0126;
Nc = 30;
c = linspace(0,Uc,Nc); 

% p = parpool(50);
for x = 1:ix
    lambdaxVector = logspace(0.5,6.5,ix);
    lambdazVector = logspace(0.5,6.5,jz);
    
    singularTmp = zeros(1,jz);
    cTmp = zeros(1,jz);
    lambdax = lambdaxVector(x);
    kx = 2*pi/lambdax*R;
    for z = 1:jz
        fprintf('x=%d,z=%d\n',x,z)
        lambdaz = lambdazVector(z);
        kz = 2*pi/lambdaz*R;
       
        singular_vector = zeros(1,Nc);
        for k = 1:Nc
            singular_vector(k) = singularvalue_up(kx,kz,c(k),N,R,flowtype,method,eddy,wall,Y);
        end
        [sinTmp,maxTmp] = max(singular_vector);
        cTmp(1,z) = c(maxTmp);
        singularTmp(1,z) = sinTmp;

    end
    singular(x,:) = singularTmp;
    cmax(x,:) =  cTmp;
end

% delete(p);
path = [method,'_',eddy,'_',wall,'/singular_N=',num2str(N),'_R=',num2str(R),'_Nc=',num2str(Nc),'_i',num2str(i)];  
save([path,'.mat'],'singular');
path = [method,'_',eddy,'_',wall,'/cmax_N=',num2str(N),'_R=',num2str(R),'_Nc=',num2str(Nc),'_i',num2str(i)]; 
save([path,'.mat'],'cmax');
