clear
clc
kx = 12;kz = 120;c = 10;  % near-wall cycle
%kx = 1;kz = 10;c = 16;   % VLSMs

flowtype = 'turbu';
method = 'IOA'; % 'IOA', 'SIOA'
wall = 'com'; % 'com', 'rigid'
eddy = 'eddyoff'; % 'eddyoff', 'eddyon'

N = 200;  % wall-normal grids
NYk = 100; % Ck grids
NY = 100; % Cd grids
R = 2000;  %Renolds number
i = sqrt(-1);
Cm = 2;
omega = c*kx;

fontn = 30;

path = [method,'_',eddy,'_',wall,'/singular_CkCd_NY=',num2str(NY),'_N=',num2str(N),'_R=',num2str(R),'_kx=',num2str(kx),'_kz=',num2str(kz),'_c=',num2str(c)];
load(path);

% colorbar scope
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff')
    str = '${\bar{\sigma}_c/\bar{\sigma}_0}$';
    if kx == 12
        cormin = 0.9;cormax = 1.1;
        ticks = [0.9,1,1.1];ticklabels={'0.9','1','1.1'};
    else
        cormin = 0.6;cormax = 1.4;
        ticks = [0.6,0.8,1,1.2,1.4];ticklabels={'0.6','0.8','1','1.2','1.4'};
    end
end
if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon')
    str = '${\bar{\sigma}_c^e /\bar{\sigma}_0^e}$';
    if kx == 12
        cormin = 0.9;cormax = 1.6;
        ticks = [0.9,1,1.2,1.4,1.6];ticklabels={'0.9','1','1.2','1.4','1.6'};
    else
        cormin = 0.8;cormax = 3;
        ticks = [0.8,1,2,3];ticklabels={'0.8','1','2.0','3.0'};
    end
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff')
    str = '${{\mu_c}/{\mu_0}}$';
    if kx == 12
        cormin = 0.95;cormax = 1.05;
        ticks = [0.95,1,1.05];ticklabels={'0.95','1','1.05'};
    else
        cormin = 0.6;cormax = 1.4;
        ticks = [0.6,0.8,1,1.2,1.4];ticklabels={'0.6','0.8','1','1.2','1.4'};
    end
end
if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon')
    str = '${{\mu_c^e}/{\mu_0^e}}$';
    if kx == 12
        cormin = 0.95;cormax = 1.05;
        ticks = [0.95,1,1.05];ticklabels={'0.95','1','1.05'};
    else
        cormin = 0.8;cormax = 3;
        ticks = [0.8,1,2,3];ticklabels={'0.8','1','2.0','3.0'};
    end
end

if kx == 12
    Ckvector = linspace(27000,31000,NYk);
    Cdvector = linspace(-10,10,NY);
end
if kx ==1
    Ckvector = linspace(350,700,NY);
    Cdvector = linspace(-10,10,NY);
end

Normrigid = singularvalue_up(kx,kz,c,N,R,flowtype,method,eddy,'rigid',[]);

Norm = Normcom./Normrigid;

pcolor(Cdvector,Ckvector,Norm);
shading interp


colormap(redblue1(cormin,cormax));
caxis([cormin,cormax]);
h=colorbar;
h.Title.String = str;
h.Title.Interpreter = 'latex';
h.Title.FontSize = fontn;
h.Ticks=ticks;
h.TickLabels=ticklabels;
hold on


set(gcf,'unit','centimeters','position',[10 7 16 14]);
set(gca,'unit','centimeters','position',[3.5,3,9,9],'fontsize',fontn,'fontname','Times')

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.015]);
set(gca,'layer','top');

xlabel('${C_d}$','Interpreter','latex');
ylabel('${C_k}$','Interpreter','latex');
grid on

Ckfun = @(Cdvar) (2.*Cm.^2.*omega.^2 + Cdvar.^2) ./ 2 ./ Cm;
Ckr = Ckfun(Cdvector);
plot(Cdvector,Ckr,'k','LineWidth',1)

[C,h] = contour(Cdvector,Ckvector,Norm,[1,1],'--k','LineWidth',1.5);
clabel(C,h,'Color','k')
hold on

% find min
[mincl,row] = min(Norm);
[minNorm,column] = min(mincl);
minCk = Ckvector(row(column));
minCd = Cdvector(column);
[ReY,ImY] = C2Y(Cm,minCk,minCd,omega);
omegar = sqrt(minCk/Cm*(1-2*(minCd/2/sqrt(minCk*Cm))^2));
fprintf('ReY=%d,ImY=%d,\nCk=%d,Cd=%d,omegar=%d,decrease=%d\n',ReY,ImY,minCk,minCd,omegar,1-minNorm);
plot(minCd,minCk,'^w','MarkerFaceColor','w','Markersize',12,'LineWidth',3)
hold on

print('-dpng','-r300', [path,'.png'])