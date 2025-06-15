clear
clc

N = 400; % wall-normal grid
R = 8900; % Reynolds number
nc = 201; % phase speed grid
nkx = 120; % kx grid
nkz = 120; % kz grid

flowtype = 'turbu'; % 'turbu'(rigid);'com3300','com6700','com8900'(com)
method = 'IOA';
eddy = 'eddyoff';
wall = 'rigid'; % 'rigid', 'com'
uvwp = 'p';
Cm = 0.46; Ck = 181; Cd=0.091;

if strcmp(wall, 'com')
    if R ==3300
        uc = 26;
    end
    if R == 6700
        uc = 28;
    end
    if R ==8900
        uc = 28;
    end
    Ck1 = Ck*uc;
    Cd1 = Cd*uc;
end

psd = zeros(nkx,nkz,nc,N);
psdmax = zeros(N,nkz,nkx);
cmax = zeros(N,nkz,nkx);

parpool(60);
parfor x = 1:nkx
    kxVector = logspace(-7,0,nkx).*R;
    kzVector = logspace(-7,0,nkz).*R;
    
    cVector = linspace(0,30,nc);
    
    psd_z = zeros(nkz,nc,N);
    psdmax_z = zeros(N,nkz);
    cmax_z = zeros(N,nkz);
    
    kx = kxVector(x);
    for z = 1:nkz
        
        kz = kzVector(z);
        psdTmp = zeros(nc,N);
        
        for cc = 1:nc
            fprintf('x=%d,z=%d,c=%d\n',x,z,cc)
            c = cVector(cc);
            
            omega = c*kx;
            if strcmp(wall, 'com')
                [ReY,ImY] = C2Y(Cm,Ck1,Cd1,omega);
                Y = ReY + sqrt(-1)*ImY;
            else
                Y=0;
            end
            
            psdTmp(cc,:) = PSD(kx,kz,c,N,R,flowtype,method,eddy,wall,Y,uvwp);
        end
        
        psd_z(z,:,:) = psdTmp;
        [psdmaxTmp,maxci] = max(psdTmp);
        cmaxTmp = cVector(maxci);
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
    path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'kplus-70_final'];
    save([path,'.mat'],'psd','-v7.3');
    path = [wall,'/',uvwp,'/psdmax_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'kplus-70_final'];
    save([path,'.mat'],'psdmax');
    path = [wall,'/',uvwp,'/cmax_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'kplus-70_final'];
    save([path,'.mat'],'cmax');
else
    path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];
    save([path,'.mat'],'psd','-v7.3');
    path = [wall,'/',uvwp,'/psdmax_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];
    save([path,'.mat'],'psdmax');
    path = [wall,'/',uvwp,'/cmax_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];
    save([path,'.mat'],'cmax');
end



