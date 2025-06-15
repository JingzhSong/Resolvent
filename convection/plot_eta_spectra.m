clear
clc

N = 400; % wall-normal grid (400 for R6700 and R8900, 300 for R3300)
R = 8900; % Reynolds number
nc = 201; % omega grid
nkx = 100;% kx grid
nkz = 100;% kz grid

method = 'IOA';
eddy = 'eddyoff';
wall = 'com';
Y = [];
Cm = 0.46; Ck = 181; Cd=0.091;
uvwp = 'v';

% %% read data and process the psd
% path = [wall,'/',uvwp,'/psd_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'Ucom_final'];
% load([path,'.mat'],'psd');
% psd = real(psd);
% psd(find(isnan(psd)==1)) = 0;
% ny = 1; % The wall displacement only needs to take the v spectrum at y=0
% 
% psd_y = squeeze(psd(:,:,:,ny));
% omegaVector = linspace(0.1,3000,nc);
% domega = omegaVector(2) - omegaVector(1);
% kxVector = linspace(0,200,nkx);
% dkx = kxVector(2)-kxVector(1);
% kzVector = linspace(0,200,nkz);
% dkz = kzVector(2)-kzVector(1);
% 
% % The v spectrum is converted into the wall displacement spectrum
% for ic = 1:nc
%     omega = omegaVector(ic);
%     psd_y(:,:,ic) = psd_y(:,:,ic) ./ omega ./ omega;
% end
% 
% kx_matrix = zeros(nkx,nkz,nc);
% kz_matrix = zeros(nkx,nkz,nc);
% for ikx = 1:nkx
%     kx_matrix(ikx,:,:) = kxVector(ikx)* ones(1,nkz,nc);
% end
% for ikz = 1:nkz
%     kz_matrix(:,ikz,:) = kzVector(ikz)* ones(nkx,1,nc);
% end
% 
% % Integrate over kz
% spec = squeeze(sum(psd_y,2)).*dkz;
% 
% PSD = (spec./(sum(spec,'all').*dkx.*domega))';
% path1 = [wall,'/',uvwp,'/processed_psd_wall_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'Ucom_final'];
% save([path1,'.mat'],'PSD');

%% plot
path1 = [wall,'/',uvwp,'/processed_psd_wall_N=',num2str(N),'nomega=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'Ucom_final'];
load([path1,'.mat'],'PSD');

fontn = 24;
[y,~] = chebdif(N,2);

if R ==3300
    uc = 26;
end
if R == 6700
    uc = 28;
end
if R ==8900
    uc = 28;
end

kxVector = linspace(0,200,nkx);
omegaVector = linspace(0.1,3000,nc);
[KX OMEGA] = meshgrid(kxVector,omegaVector);

figure;
pcolor(KX,OMEGA,log10(PSD)); % 整体积分是1，归一化
% pcolor(KX,OMEGA,log10(spec)'); % 不归一化
xlabel('$k_x$','Interpreter','latex');
ylabel('$\omega$','Interpreter','latex');
hold on

shading interp
load('MyColormaps.mat')
colormap(mycmap);
caxis([-12,-2])
h = colorbar('Location', 'northoutside');
h.Label.String = '${\rm log}_{10} \Phi_{\hat{\eta}}$'; 
h.Label.Interpreter = 'latex';                      
h.Label.FontSize = fontn;             

axis([0 200 0 3000])

% Find the maximum of different omega values
[maxspec,nmax] = max(PSD(1:50,:));
kxx = kxVector;
hold on

% Plot the wall convective velocity
plot(kxVector,0.44 .*uc.* kxVector,'.k','LineWidth',2)
hold on
plot(kxVector,0.53 .*uc.* kxVector,'--k','LineWidth',2)
hold on

% Plot resonance frequency
yline(19.78.*uc,'k','LineWidth',2)

set(gca,'layer','top')

h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;

set(gcf,'unit','centimeters','position',[10 7 23 15]);
set(gca,'unit','centimeters','position',[3,2.5,13,9],'fontsize',fontn,'fontname','Times')

% dual y axis
ax1 = gca;
ax2 = axes('Position', ax1.Position, ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...
           'XTick', [], ...
           'Box', 'off'); 
ylim_right = ax1.YLim / uc;
ax2.YLim = ylim_right;
ylabel(ax2, '$\omega / U_c^+$', 'Interpreter', 'latex', 'FontSize', fontn);
ax2.YTick = 0:20:100;

set(gcf,'unit','centimeters','position',[10 7 20 15]);
set(gca,'unit','centimeters','position',[3,2.5,13,9],'fontsize',fontn,'fontname','Times')

print('-dpng','-r300', [path1,'.png'])