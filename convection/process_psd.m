clear
clc

N = 400; % wall-normal grid
R = 8900; % Reynolds number
nc = 201; % omega grid
nkx = 120; % kx grid
nkz = 120; % kz grid

method = 'IOA';
eddy = 'eddyoff';
wall = 'com';
uvwp = 'p';
Y = [];

% Small scales are not considered
lambdax_min = 500;
lambdaz_min = 80;

kxVector = logspace(-7,0,nkx).*R;
kzVector = logspace(-7,0,nkz).*R;

lambdaxVector = 2.*pi./kxVector .*R;
lambdazVector = 2.*pi./kzVector .*R;

if strcmp(wall, 'rigid')
    path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];
else
    Cm = 0.46; Ck = 181; Cd=0.091;
    path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'kplus-70_final'];
end

%% read data
path0 = [wall,'/',uvwp];

path = [path0,'/psd',path1];
load([path,'.mat'],'psd');
path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
psd = real(psd);

psd(find(isnan(psd))) = 0;
psd_ave = squeeze(sum(psd,3)); % sum of omega
%%
path = [path0,'/psd_kx_kz',path1];
save([path,'.mat'],'psd_ave');