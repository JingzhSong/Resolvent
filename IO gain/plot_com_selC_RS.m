clear
clc
kx = 12;kz = 120;c = 10;  % near-wall cycle
%kx = 1;kz = 10;c = 16;   % VLSMs

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
wall = 'com'; % 'com', 'rigid'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'

N = 200;  % wall-normal grids
NY = 100; % Ck and Cd grids
R = 2000;  %Renolds number
i = sqrt(-1);
Cm = 2;
omega = c*kx;
[y,~] = chebdif(N,2);

fontn = 30;

path = [method,'_',eddy,'_',wall,'/RS_CkCd_NY=',num2str(NY),'_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kx),'_kz=',num2str(kz),'_c=',num2str(c)];
load(path);

if kx == 12
Ckvector = linspace(27000,31000,NY);
Cdvector = linspace(-10,10,NY);
end
if kx == 1
Ckvector = linspace(350,700,NY);
Cdvector = linspace(-10,10,NY);
end

[~,~,~,~,~,~,~,~,RSrigid] = UF_up(kx,kz,c,N,R,flowtype,method,eddy,'rigid',[]);

% colorbar scope
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff')
    if kx == 12
        cormin = 0.8;cormax = 1.2;
        ticks = [0.8,0.9,1,1.1,1.2];ticklabels={'0.8','0.9','1','1.1','1.2'};
    else
        cormin = 0;cormax = 2;
        ticks = [0,1,2];ticklabels={'0','1','2'};
    end
end
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon')
    if kx == 12
        cormin = 0.8;cormax = 1.2;
        ticks = [0.8,0.9,1,1.1,1.2];ticklabels={'0.8','0.9','1','1.1','1.2'};
    else
        cormin = 0;cormax = 2;
        ticks = [0,1,2];ticklabels={'0','1','2'};
    end
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff')
    if kx == 12
        cormin = -4;cormax = 6;
        ticks = [-4,1,6];ticklabels={'-4','1','6'};
    else
        cormin = -100;cormax = 50;
        ticks = [-100,1,50];ticklabels={'-100','1','50'};
    end
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon')
    if kx == 12
        cormin = -4;cormax = 6;
        ticks = [-4,1,6];ticklabels={'-4','1','6'};
    else
        cormin = -100;cormax = 50;
        ticks = [-100,1,50];ticklabels={'-100','1','50'};
    end
end


RS = RScom./RSrigid;

pcolor(Cdvector,Ckvector,RS);
shading interp
caxis([cormin,cormax]);

colormap(redblue1(cormin,cormax));
caxis([cormin,cormax]);
h=colorbar;
h.Title.String = '${RS_c/RS_0}$';
h.Title.Interpreter = 'latex';
h.Title.FontSize = fontn;
h.Ticks=ticks;
h.TickLabels=ticklabels;
hold on

set(gcf,'unit','centimeters','position',[10 7 16 14]);
set(gca,'unit','centimeters','position',[3.5,3,9,9],'fontsize',fontn,'fontname','Times')

set(gca,'xminortick','on');
set(gca,'ticklength',[0.025 0.01]);
set(gca,'Layer','top');
grid off;

xlabel('${C_d}$','Interpreter','latex');
ylabel('${C_k}$','Interpreter','latex');

Ckfun = @(Cdvar) (2.*Cm.^2.*omega.^2 + Cdvar.^2) ./ 2 ./ Cm;
Ckr = Ckfun(Cdvector);
plot(Cdvector,Ckr,'k','LineWidth',1)

[C,h] = contour(Cdvector,Ckvector,RS,[1,1],'--k','LineWidth',1.5);
clabel(C,h,'Color','k')
hold on

% find min
[mincl,row] = min(RS);
[minNorm,column] = min(mincl);
minCk = Ckvector(row(column));
minCd = Cdvector(column);
[ReY,ImY] = C2Y(Cm,minCk,minCd,omega);
omegar = sqrt(minCk/Cm*(1-2*(minCd/2/sqrt(minCk*Cm))^2));
fprintf('ReY=%d,ImY=%d,\nCk=%d,Cd=%d,\nomegar=%d,decrease=%d\n',ReY,ImY,minCk,minCd,omegar,1-minNorm);
plot(minCd,minCk,'^w','MarkerFaceColor','w','Markersize',12,'LineWidth',3)
hold on

print('-dpng','-r300', [path,'.png'])