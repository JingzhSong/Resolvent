clear
clc

kz = 0;
N = 200;  % wall-normal grids
R = 140;  %Renolds number
k = 0.426; alpha = 25.4;
i = sqrt(-1);

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'
wall = 'rigid';

Cm = 2; Ck = 440; Cd = 10.5;
omegar = 14.36;

NY = 50;
kxvector = logspace(0,1,NY);
cvector = linspace(4,14,NY);
NormRigid = zeros(NY,NY);


%% 计算不同波数
% p = parpool(50);
for cc = 1:NY
    Normtmp = zeros(1,NY);
    c = cvector(cc);
    cc
    for kxx = 1:NY
        kx = kxvector(kxx);
        kxx
        Y = 0;
        
        Normtmp(kxx) = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);
    end
    NormRigid(cc,:) = Normtmp;
end
% delete(p);
path = [method,'_',eddy,'_',wall,'/Norm_DNS_2D'];
save([path,'.mat'],'NormRigid');
