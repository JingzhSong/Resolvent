%% IOA fix kz
clear
clc
N=200;  % wall-normal grids
R = 2000;  %Renolds number
i=50;  % kx grids
j=50;  % kz grids
fontn = 35;
lambdaxVector = logspace(0.5,6.5,i);
lambdazVector= logspace(0.5,6.5,j);

method = 'IOA';
wall = 'rigid';
eddy = 'eddyoff';
lambdaz = 100;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl1 = singular(:,1);
[sgl1max,n_sgl1max] = max(sgl1);

lambdaz = 0.6*R;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl2 = singular(:,1);
[sgl2max,n_sgl2max] = max(sgl2);

eddy = 'eddyon';
lambdaz = 100;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl3 = singular(:,1);
[sgl3max,n_sgl3max] = max(sgl3);

lambdaz = 0.6*R;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl4 = singular(:,1);
[sgl4max,n_sgl4max] = max(sgl4);

figure;

loglog(lambdaxVector,sgl1,'b','LineWidth',4);
hold on
loglog(lambdaxVector,sgl2,'--k','LineWidth',4);
hold on
loglog(lambdaxVector,sgl3,'ob','LineWidth',2);
hold on
loglog(lambdaxVector,sgl4,'sk','LineWidth',2); 

loglog([lambdaxVector(n_sgl1max),lambdaxVector(n_sgl1max)],[10^(-0.3),sgl1max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl3max),lambdaxVector(n_sgl3max)],[10^(-0.3),sgl3max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl2max),lambdaxVector(n_sgl2max)],[sgl2max,10^3.5],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl4max),lambdaxVector(n_sgl4max)],[sgl4max,10^3.5],'--k','LineWidth',1.5)

xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('$\Vert \mathcal{H} \Vert_{\infty}, \Vert \mathcal{H}^{e} \Vert_{\infty}$','Interpreter','latex');
legend('$\Vert \mathcal{H} \Vert_{\infty}, \lambda_z^+=100$','$\Vert \mathcal{H} \Vert_{\infty}, \lambda_z=0.6$',...
    '$\Vert \mathcal{H}^{e} \Vert_{\infty},\lambda_z^+=100$','$\Vert \mathcal{H}^{e} \Vert_{\infty},\lambda_z=0.6$',...
    'location','NorthWest','Interpreter','latex','fontsize',fontn);

axis([10^0 10^6 10^(-4) 10^5])
xticks(logspace(0,6,7))
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')
box off;
hold on

% dual x axis
ax1 = gca; 
ax1.LineWidth = 1; 
ax1 = axes('Position', ax1.Position, 'TickDir','in','XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2 = axes('Position', ax1.Position, 'TickDir','out','XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2.XLim = [10^0/R 10^6/R];
ax2.YLim = [10^(-4) 10^5];
ax2.YTickLabel = [];
xticks(logspace(-3,3,7))
xlabel('${\lambda}_x$','Interpreter','latex');
ax2.LineWidth = 1;
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')

path = [method,'_','eddyon','_',wall,'/singular_line_lambdaz_N=',num2str(N),'_R=',num2str(R),'new'];

print('-dpng','-r300', [path,'.png'])



%% SIOA fix kz
clear
clc
N = 200;
R = 2000;
i = 50;
j = 50;
fontn = 35;
lambdaxVector = logspace(0.5,6.5,i);

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyoff';
wall = 'rigid';
lambdaz = 100;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl1 = singular(:,1);
[sgl1max,n_sgl1max] = max(sgl1);

lambdaz = 0.6*R;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl2 = singular(:,1);
[sgl2max,n_sgl2max] = max(sgl2);

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyon';
wall = 'rigid';
lambdaz = 100;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl3 = singular(:,1);
[sgl3max,n_sgl3max] = max(sgl3);

lambdaz = 0.6*R;
path = [method,'_',eddy,'_',wall,'/singular_lambdaz=',num2str(lambdaz),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl4 = singular(:,1);
[sgl4max,n_sgl4max] = max(sgl4);

figure;

loglog(lambdaxVector,sgl1,'b','LineWidth',3);
hold on
loglog(lambdaxVector,sgl2,'--k','LineWidth',3);
hold on
loglog(lambdaxVector,sgl3,'ob','LineWidth',1.5);
hold on
loglog(lambdaxVector,sgl4,'sk','LineWidth',1.5);
hold on

lambdaxVector(n_sgl1max)
lambdaxVector(n_sgl2max)
lambdaxVector(n_sgl3max)
lambdaxVector(n_sgl4max)
loglog([lambdaxVector(n_sgl1max),lambdaxVector(n_sgl1max)],[10^(-0.3),sgl1max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl3max),lambdaxVector(n_sgl3max)],[10^(-0.3),sgl3max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl2max),lambdaxVector(n_sgl2max)],[sgl2max,10^3.5],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl4max),lambdaxVector(n_sgl4max)],[sgl4max,10^3.5],'--k','LineWidth',1.5)


xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('$\Vert \mathcal{H}_\nabla\Vert_{\mu}, \Vert \mathcal{H}^{e}_\nabla \Vert_{\mu}$','Interpreter','latex');
legend('$\Vert \mathcal{H}_\nabla \Vert_{\mu}, \lambda_z^+=100$','$\Vert \mathcal{H}_\nabla \Vert_{\mu}, \lambda_z=0.6$',...
    '$\Vert \mathcal{H}^{e}_\nabla \Vert_{\mu},\lambda_z^+=100$','$\Vert \mathcal{H}^{e}_\nabla \Vert_{\mu},\lambda_z=0.6$',...
'location','NorthWest','Interpreter','latex','fontsize',fontn);

axis([10^0 10^6 10^(-0.3) 10^3.5])
xticks(logspace(0,6,7))
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')

box off;
hold on

% dual x axis
ax1 = gca;
ax1.LineWidth = 1;
ax1 = axes('Position', ax1.Position, 'TickDir','in','XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2 = axes('Position', ax1.Position, 'TickDir','out','XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2.XLim = [10^0/R 10^6/R];
ax2.YLim = [10^(-0.3) 10^3.5];
ax2.YTickLabel = [];
xticks(logspace(-3,3,7))
xlabel('${\lambda}_x$','Interpreter','latex');
ax2.LineWidth = 1;
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')

path = [method,'_','eddyon','_',wall,'/singular_line_lambdaz_N=',num2str(N),'_R=',num2str(R),'new'];

print('-dpng','-r300', [path,'.png'])




%% IOA fix kx
clear
clc
N=200;
R = 2000;
i=50;
j=50;
fontn = 35;
lambdaxVector = logspace(0.5,6.5,i);
lambdazVector= logspace(0.5,6.5,j);

method = 'IOA';
wall = 'rigid';
eddy = 'eddyoff';
lambdax = 1000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl1 = singular(:,1);
[sgl1max,n_sgl1max] = max(sgl1);

lambdax = 12000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl3 = singular(:,1);
[sgl3max,n_sgl3max] = max(sgl3);

eddy = 'eddyon';
lambdax = 1000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl4 = singular(:,1);
[sgl4max,n_sgl4max] = max(sgl4);

lambdax = 12000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl6 = singular(:,1);
[sgl6max,n_sgl6max] = max(sgl6);

figure;

loglog(lambdaxVector,sgl1,'b','LineWidth',3);
hold on
loglog(lambdaxVector,sgl3,'--k','LineWidth',3);
hold on
loglog(lambdaxVector,sgl4,'ob','LineWidth',1.5);
hold on
loglog(lambdaxVector,sgl6,'sk','LineWidth',1.5);


lambdaxVector(n_sgl1max)
lambdaxVector(n_sgl3max)
lambdaxVector(n_sgl4max)
lambdaxVector(n_sgl6max)

loglog([lambdaxVector(n_sgl1max),lambdaxVector(n_sgl1max)],[10^(-4),sgl1max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl4max),lambdaxVector(n_sgl4max)],[10^(-4),sgl4max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl3max),lambdaxVector(n_sgl3max)],[sgl3max,10^6],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl6max),lambdaxVector(n_sgl6max)],[sgl6max,10^6],'--k','LineWidth',1.5)

xlabel('${\lambda}_z^+$','Interpreter','latex');
ylabel('$\Vert \mathcal{H} \Vert_{\infty}, \Vert \mathcal{H}^{e} \Vert_{\infty}$','Interpreter','latex');
legend('$\Vert \mathcal{H} \Vert_{\infty}, \lambda_x^+=1000$','$\Vert \mathcal{H} \Vert_{\infty}, \lambda_x=6$',...
    '$\Vert \mathcal{H}^{e} \Vert_{\infty},\lambda_x^+=1000$','$\Vert \mathcal{H}^{e} \Vert_{\infty},\lambda_x=6$',...
    'location','NorthWest','Interpreter','latex','fontsize',fontn);
axis([10^0 10^6 10^(-4) 10^6])
xticks(logspace(0,6,7))

set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')
box off;
hold on

% dual x axis
ax1 = gca;
ax1.LineWidth = 1;
ax1 = axes('Position', ax1.Position, 'TickDir','in','XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2 = axes('Position', ax1.Position, 'TickDir','out','XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2.XLim = [10^0/R 10^6/R]; 
ax2.YLim = [10^(-4) 10^6];
ax2.YTickLabel = [];
xticks(logspace(-3,3,7))
xlabel('${\lambda}_z$','Interpreter','latex');
ax2.LineWidth = 1;
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')

path = [method,'_','eddyon','_',wall,'/singular_line_lambdax_N=',num2str(N),'_R=',num2str(R)];

print('-dpng','-r300', [path,'new.png'])



%% SIOA fix kx
clear
clc
N = 200;
R = 2000;
i = 50;
j = 50;
fontn = 35;
lambdaxVector = logspace(0.5,6.5,i);

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyoff';
wall = 'rigid';
lambdax = 1000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl1 = singular(:,1);
[sgl1max,n_sgl1max] = max(sgl1);

lambdax = 12000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl3 = singular(:,1);
[sgl3max,n_sgl3max] = max(sgl3);

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyon';
wall = 'rigid';
lambdax = 1000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl4 = singular(:,1);
[sgl4max,n_sgl4max] = max(sgl4);

lambdax = 12000;
path = [method,'_',eddy,'_',wall,'/singular_lambdax=',num2str(lambdax),'_N=',num2str(N),'_R=',num2str(R)];
load([path,'.mat'],'singular');
sgl6 = singular(:,1);
[sgl6max,n_sgl6max] = max(sgl6);

figure;

loglog(lambdaxVector,sgl1,'b','LineWidth',3);
hold on
loglog(lambdaxVector,sgl3,'--k','LineWidth',3);
hold on
loglog(lambdaxVector,sgl4,'ob','LineWidth',1.5);
hold on
loglog(lambdaxVector,sgl6,'sk','LineWidth',1.5);
hold on

lambdaxVector(n_sgl1max)
lambdaxVector(n_sgl3max)
lambdaxVector(n_sgl4max)
lambdaxVector(n_sgl6max)

loglog([lambdaxVector(n_sgl1max),lambdaxVector(n_sgl1max)],[10^(-0.3),sgl1max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl4max),lambdaxVector(n_sgl4max)],[10^(-0.3),sgl4max],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl3max),lambdaxVector(n_sgl3max)],[sgl3max,10^4],'--k','LineWidth',1.5)
loglog([lambdaxVector(n_sgl6max),lambdaxVector(n_sgl6max)],[sgl6max,10^4],'--k','LineWidth',1.5)

xlabel('${\lambda}_z^+$','Interpreter','latex');
ylabel('$\Vert \mathcal{H}_\nabla\Vert_{\mu}, \Vert \mathcal{H}^{e}_\nabla \Vert_{\mu}$','Interpreter','latex');
legend('$\Vert \mathcal{H}_\nabla \Vert_{\mu}, \lambda_x^+=1000$','$\Vert \mathcal{H}_\nabla \Vert_{\mu}, \lambda_x=6$',...
    '$\Vert \mathcal{H}^{e}_\nabla \Vert_{\mu},\lambda_x^+=1000$','$\Vert \mathcal{H}^{e}_\nabla \Vert_{\mu},\lambda_x=6$',...
'location','NorthWest','Interpreter','latex','fontsize',fontn);
axis([10^0 10^6 10^(-0.3) 10^4])
xticks(logspace(0,6,7))
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')
box off;

% dual x axis
ax1 = gca;
ax1.LineWidth = 1;
ax1 = axes('Position', ax1.Position, 'TickDir','in','XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2 = axes('Position', ax1.Position, 'TickDir','out','XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XScale', 'log', 'YScale', 'log');
ax2.XLim = [10^0/R 10^6/R];
ax2.YLim = [10^(-0.3) 10^4];
ax2.YTickLabel = [];
xticks(logspace(-3,3,7))
xlabel('${\lambda}_z$','Interpreter','latex');
ax2.LineWidth = 1;
set(gcf,'unit','centimeters','position',[8 2 32 25]);
set(gca,'unit','centimeters','position',[4.5,4.5,25,16],'fontsize',fontn,'fontname','Times')

path = [method,'_','eddyon','_',wall,'/singular_line_lambdax_N=',num2str(N),'_R=',num2str(R)];

print('-dpng','-r300', [path,'.png'])
