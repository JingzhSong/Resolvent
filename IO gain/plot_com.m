%% 绘图
clear
clc
N = 200; % wall-normal grids
R = 2000; %Renolds number
i=50; % kx grid
j=50; %kz grid

Nc = 30;
kx = 12;kz = 120;c = 10;  % near-wall cycle
%kx = 1;kz = 10;c = 16;   % VLSMs

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
wall = 'com';
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'
selc = 'SIN'; % 'SIN', 'RS' - objective functions
effect = 'good';

lambdaxVector = logspace(0.5,6.5,i);
lambdazVector= logspace(0.5,6.5,j);

% read data for compliant
path = [method,'_',eddy,'_',wall,'/singular_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kx),'_ikx=',num2str(j),'_Nc=',num2str(Nc),selc,'_',effect];  
load([path,'.mat']);
path = [method,'_',eddy,'_',wall,'/cmax_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kx),'_ikx=',num2str(j),'_Nc=',num2str(Nc),selc,'_',effect];  
load([path,'.mat']);

singular_com = singular;

% read data for rigid
path = [method,'_',eddy,'_rigid/singular_N=',num2str(N),'_R=',num2str(R),'_Nc=',num2str(Nc),'_i',num2str(i)];
load([path,'.mat']);
singular_rigid = singular;

%% plot compliant wall
figure;
fontn = 20;

if strcmp(method,'SIOA')
    cormin = 0;
    if strcmp(eddy,'eddyoff')
        cormax = 4;
    else
        cormax = 3;
    end
else
    cormin = -4;
    if strcmp(eddy,'eddyoff')
        cormax = 6;
    else
        cormax = 2;
    end
end

A = log10(singular_com');

pcolor(lambdazVector,lambdaxVector,A);
set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');

shading interp

if strcmp(method, 'IOA') && strcmp(eddy,'eddyoff')
     colorbar_title = '${\rm log_{10}}\Vert \mathcal{H} \Vert_{\infty}$';
end
if strcmp(method, 'SIOA') && strcmp(eddy,'eddyoff')
     colorbar_title = '${\rm log_{10}}\Vert \mathcal{H}_{\nabla} \Vert_{\mu}$';
end
if strcmp(method, 'IOA') && strcmp(eddy,'eddyon')
     colorbar_title = '${\rm log_{10}}\Vert \mathcal{H}^{e} \Vert_{\infty}$';
end
if strcmp(method, 'SIOA') && strcmp(eddy,'eddyon')
     colorbar_title = '${\rm log_{10}}\Vert \mathcal{H}^{e}_{\nabla} \Vert_{\mu}$';
end

load('MyColormaps.mat')
colormap(mycmap);

h=colorbar;
title(colorbar_title,'Interpreter','latex');
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;
if strcmp(method, 'IOA')
    h.Ticks=[-4,-2,0,2,4,6];
    h.TickLabels={'-4','-2','0','2','4','6'};
else
    h.Ticks=[0,1,2,3,4];
    h.TickLabels={'0','1','2','3','4'};
end

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')

caxis([cormin,cormax]);
set(gcf,'unit','centimeters','position',[10 7 14 13.5]);
set(gca,'unit','centimeters','position',[3,3,9,9],'fontsize',fontn,'fontname','Times')
hold on
plot(2*pi/kx*R,2*pi/kz*R,'ko','MarkerSize',10,'LineWidth',3)

path = [method,'_',eddy,'_',wall,'/singular_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kx),'_ikx=',num2str(j),'_Nc=',num2str(Nc),selc,'_',effect];  
saveas(gcf, [path,'.png']);

