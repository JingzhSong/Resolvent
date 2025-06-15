clear
clc

N = 200;  % wall-normal grids
R = 2000; %Renolds number
ix=50; % kx grid
jz=50; % kz grid
nt=ix*jz;
Uc= 24.0126; % centerline velocity
Nc = 30; % c grid

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
wall = 'com'; % 'com', 'rigid'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'
selc = 'SIN'; % 'SIN', 'RS' - objective function
effect = 'good';

% kxsel = 12; %near-wall cycle
kxsel = 1; %VLSMs

i = sqrt(-1);

singular = zeros(ix,jz);
cmax = zeros(ix,jz);

%%
% p = parpool(50);
for x = 1:ix
    lambdaxVector = logspace(0.5,6.5,ix);
    lambdazVector = logspace(0.5,6.5,jz);
    cvecter = linspace(0,Uc,Nc); 
    
    singularTmp = zeros(1,jz);
    cTmp = zeros(1,jz);
    lambdax = lambdaxVector(x);
    kx = 2*pi/lambdax*R;
    for z = 1:jz
        fprintf('x=%d,z=%d\n',x,z)
        lambdaz = lambdazVector(z);
        kz = 2*pi/lambdaz*R;
       
        singular_vector = zeros(1,Nc);
        for k = 2:Nc
            c = cvecter(k);
            omega=c*kx;
            [Cm,Cd,Ck,Y] = readcom_final(method,eddy,wall,selc,effect,kxsel,omega);
            singular_vector(k) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
        end
        [sinTmp,maxTmp] = max(singular_vector);
        cTmp(1,z) = cvecter(maxTmp);
        singularTmp(1,z) = sinTmp;

    end
    singular(x,:) = singularTmp;
    cmax(x,:) =  cTmp;
end

% delete(p);
path = [method,'_',eddy,'_',wall,'/singular_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kxsel),'_ikx=',num2str(ix),'_Nc=',num2str(Nc),selc,'_',effect];  
save([path,'.mat'],'singular');
path = [method,'_',eddy,'_',wall,'/cmax_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kxsel),'_ikx=',num2str(ix),'_Nc=',num2str(Nc),selc,'_',effect]; 
save([path,'.mat'],'cmax');
