clear
clc

N = 400; % wall-normal grid
R = 8900; % Reynolds number
nc = 201; % omega grid
nkx = 100; % kx grid
nkz = 100; % kz grid

flowtype = 'com8900'; % 'com3300','com6700','com8900'(for different Reynolds number)
method = 'IOA';
eddy = 'eddyoff';
wall = 'com';

uvwp = 'v';

if strcmp(wall, 'com')
    Cm = 0.46; Ck = 181; Cd=0.091;
if R ==3300
    uc = 26;
end
if R == 6700
    uc = 28;
end
if R ==8900
    uc = 28;
end  
    Ck1 = Ck*uc^2;
    Cd1 = Cd*uc;
end

omegaVector = linspace(0.1,3000,nc);

psd = zeros(nkx,nkz,nc,N);
psdmax = zeros(N,nkz,nkx);
cmax = zeros(N,nkz,nkx);

parpool(60);
parfor x = 1:nkx
    kxVector = linspace(0,200,nkx);
    kzVector = linspace(0,200,nkz);
    psd_z = zeros(nkz,nc,N);
    psdmax_z = zeros(N,nkz);
    cmax_z = zeros(N,nkz);
    
    kx = kxVector(x);
    for z = 1:nkz
        
        kz = kzVector(z);
        psdTmp = zeros(nc,N);
        
        for cc = 1:nc
            fprintf('x=%d,z=%d,c=%d\n',x,z,cc)
            omega = omegaVector(cc);

            c = omega / kx;
            if strcmp(wall, 'com')
            [ReY,ImY] = C2Y(Cm,Ck1,Cd1,omega); % Y是omega的函数
            Y = ReY + sqrt(-1)*ImY;
            else
                Y=0;
            end
            psdTmp(cc,:) = PSD(kx,kz,c,N,R,flowtype,method,eddy,wall,Y,uvwp);
        end
        
        psd_z(z,:,:) = psdTmp;
        [psdmaxTmp,maxci] = max(psdTmp);
        cmaxTmp = omegaVector(maxci)/kx;
        psdmax_z(:,z) = psdmaxTmp';
        cmax_z(:,z) = cmaxTmp';
    end
    psd(x,:,:,:) = psd_z;
    psdmax(:,:,x) = psdmax_z;
    cmax(:,:,x) = cmax_z;
end
delete(gcp('nocreate'));

psdmax = permute(psdmax,[3,2,1]);
cmax = permute(cmax,[3,2,1]);

%%
if strcmp(wall, 'com')
    path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'Ucom_final'];
    save([path,'.mat'],'psd','-v7.3');
    path = [wall,'/',uvwp,'/psdmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'Ucom_final'];
    save([path,'.mat'],'psdmax');
    path = [wall,'/',uvwp,'/cmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'Ucom_final'];
    save([path,'.mat'],'cmax');
else
    path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'Ucom_final'];
    save([path,'.mat'],'psd','-v7.3');
    path = [wall,'/',uvwp,'/psdmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'Ucom_final'];
    save([path,'.mat'],'psdmax');
    path = [wall,'/',uvwp,'/cmax_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'Ucom_final'];
    save([path,'.mat'],'cmax');
end

