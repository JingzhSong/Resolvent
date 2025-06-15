clear
clc

N = 200; % wall-normal grid
R = 2000; % Reynolds number
nc = 201; % omega grid
nkx = 120; % kx grid
nkz = 120; % kz grid

method = 'IOA';
eddy = 'eddyoff';
wall = 'rigid';
Y = [];

% Small scales are not considered
lambdax_min = 500;
lambdaz_min = 80;

kxVector = logspace(-7,0,nkx).*R;
kzVector = logspace(-7,0,nkz).*R;

lambdaxVector = 2.*pi./kxVector .*R;
lambdazVector = 2.*pi./kzVector .*R;
kx_matrix = zeros(nkx,nkz,N);
kz_matrix = zeros(nkx,nkz,N);
for ikx = 1:nkx
    kx_matrix(ikx,:,:) = kxVector(ikx)* ones(1,nkz,N);
end
for ikz = 1:nkz
    kz_matrix(:,ikz,:) = kzVector(ikz)* ones(nkx,1,N);
end
n1=length(find(lambdaxVector > lambdax_min));
fprintf('lambdax = %d \n',lambdaxVector(n1))
n2=length(find(lambdazVector > lambdaz_min));
fprintf('lambdaz = %d \n',lambdazVector(n2))
path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];


%% read data
uvwp = 'u';
path0 = [wall,'/',uvwp];
path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');

h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix; % The weighted log integral is kx times kz. And d\omega=dc*kx, should times another kx.
u_vec = squeeze(sum(cmax(1:n1,1:n2,N/2:end) .* h(1:n1,1:n2,N/2:end),[1 2])) ./ squeeze(sum(h(1:n1,1:n2,N/2:end),[1 2])); % Small scales are not considered

%% read data
uvwp = 'v';
path0 = [wall,'/',uvwp];
path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');

h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix; % The weighted log integral is kx times kz. And d\omega=dc*kx, should times another kx.
v_vec = squeeze(sum(cmax(1:n1,1:n2,N/2:end) .* h(1:n1,1:n2,N/2:end),[1 2])) ./ squeeze(sum(h(1:n1,1:n2,N/2:end),[1 2])); % Small scales are not considered

%% read data
uvwp = 'w';
path0 = [wall,'/',uvwp];
path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');

h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix; 
w_vec = squeeze(sum(cmax(1:n1,1:n2,N/2:end) .* h(1:n1,1:n2,N/2:end),[1 2])) ./ squeeze(sum(h(1:n1,1:n2,N/2:end),[1 2]));

%%
uvwp = 'p';
path0 = [wall,'/',uvwp];
path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');
path = [path0,'/psd_kx_kz',path1];
load([path,'.mat'],'psd_ave');

h = psd_ave.*kx_matrix.*kx_matrix.*kx_matrix.*kx_matrix.*kz_matrix;
p_vec = squeeze(sum(cmax(1:n1,1:n2,N/2:end) .* h(1:n1,1:n2,N/2:end),[1 2])) ./ squeeze(sum(h(1:n1,1:n2,N/2:end),[1 2]));

%% u mean
[y,~] = chebdif(N,2);
k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
DUDy = @(y) R.*(-y)./NuT(y);
U0 = zeros(N,1);
for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end


%% plot
figure;
nline = 3;
[y,~] = chebdif(N,2);
yplus = (y+1)*R;
n=10; % Close to the center of the channel is not drawn
semilogx(yplus(N/2+n:end,1),u_vec(1+n:end),'--b','linewidth',nline); hold on
semilogx(yplus(N/2+n:end,1),v_vec(1+n:end),'--','Color',[0.196 0.804 0.196],'linewidth',nline); hold on
semilogx(yplus(N/2+n:end,1),w_vec(1+n:end),'--m','linewidth',nline); hold on
semilogx(yplus(N/2+n:end,1),p_vec(1+n:end),'r','linewidth',nline); hold on

% % DNS
nsize = 14;
load('DNSp.mat','DNSp')
load('DNSu.mat','DNSu')
load('DNSv.mat','DNSv')
load('DNSw.mat','DNSw')
semilogx(DNSu(1:2:end,1),DNSu(1:2:end,2),'bo','linewidth',nline,'Markersize',nsize);hold on
semilogx(DNSv(1:2:end,1),DNSv(1:2:end,2),'*','Color',[0.196 0.804 0.196],'linewidth',nline,'Markersize',nsize);hold on
semilogx(DNSw(1:2:end,1),DNSw(1:2:end,2),'mx','linewidth',nline,'Markersize',nsize);hold on
semilogx(DNSp(1:2:end,1),DNSp(1:2:end,2),'r^','linewidth',nline,'Markersize',nsize);hold on

semilogx(yplus(N/2:end,1),U0(N/2:end,1),':k','linewidth',nline);hold on
grid on
fontn = 28;
axis([10^(0) 10^3.5 0 25])

xlabel('$y^+$','Interpreter','latex');
ylabel('$\bar{c}^+_\psi,U^+$','Interpreter','latex');

set(gcf,'unit','centimeters','position',[10 7 26 17]);
set(gca,'unit','centimeters','position',[3,3,22,13.5],'fontsize',fontn,'fontname','Times')
legend('$\bar{c}^+_u$','$\bar{c}^+_v$','$\bar{c}^+_w$','$\bar{c}^+_p$','DNS - $\bar{c}^+_u$','DNS - $\bar{c}^+_v$','DNS - $\bar{c}^+_w$','DNS - $\bar{c}^+_p$',...
   '$U^+$','Interpreter','latex','location','SouthEast')

path = [path0,'/cov',path1,num2str(lambdax_min)];
print('-dpng','-r300', [path,'.png'])