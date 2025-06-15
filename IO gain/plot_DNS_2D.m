clear
clc
kz = 0;
N = 200; % wall-normal grids
R = 140; %Renolds number

Cm = 2; Ck = 440; Cd = 10.5;
omegar = 14.36;

flowtype = 'turbu'; 
method = 'IOA'; % 'IOA', 'SIOA'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'

NY = 50;
kxvector = logspace(0,1,NY);
cvector = linspace(4,14,NY);
NormCom = zeros(NY,NY);
    
wall = 'com';
path = [method,'_',eddy,'_',wall,'/Norm_DNS_2D.mat'];
load(path)
Norm_com = NormCom;
wall = 'rigid';
path = [method,'_',eddy,'_',wall,'/Norm_DNS_2D.mat'];
load(path)
Norm_rigid = NormCom;

% plot
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\bar{\sigma}_c}/{\bar{\sigma}_0})$';
end
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\bar{\sigma}_c^e}/{\bar{\sigma}_0^e})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff')
    str = '${\rm log_{10}}({\mu_c}/{\mu_0})$';
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon')
    str = '${\rm log_{10}}({\mu_c^e}/{\mu_0^e})$';
end

cormin=0;
cormax=1.2;

contourf(kxvector,cvector,log10(Norm_com./Norm_rigid),[0:0.1:cormax],'EdgeColor','none');
colorbar
caxis([cormin,cormax]);

set(gca,'position',[0.21,0.17,0.32*1.7,0.448*1.7],'fontsize',18,'fontname','Times');
set(gca,'XScale','log');
xlabel('$k_{x}$','Interpreter','latex');
ylabel('$c^+$','Interpreter','latex');
set(gca,'XTick',10.^(0:0.5:1),'XTickLabel',{'$10^{0}$','$10^{0.5}$','$10^{1}$'},'TickLabelInterpreter','latex');
set(gca,'YTick',[4:2:14],'YTickLabel',{'$4$','$6$','$8$','$10$','$12$','$14$'},'TickLabelInterpreter','latex');
axis([1,10,4,14]);

load('MyColormaps.mat')
colormap(mycmap);

h=colorbar;
fontn =20;
title(str,'Interpreter','latex');
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;
h.Ticks=0:0.2:cormax;

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
set(gcf,'unit','centimeters','position',[10 7 14 12]);
set(gca,'unit','centimeters','position',[2,2,9,9],'fontsize',fontn,'fontname','Times')
hold on

x = kxvector;
y = omegar./x;
plot(kxvector,y,'k','LineWidth',2);
hold on

plot(2.6,4.78,'ok','MarkerSize',10,'LineWidth',3);

wall = 'com';
path = [method,'_',eddy,'_',wall,'/Norm_DNS_2D'];
saveas(gcf, [path,'.png']);