clear
clc
N = 200;
R = 2000;
i=50;
j=50;
Nc = 30;

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyoff';
wall = 'rigid';

fontn = 20;
lambdaxVector = logspace(0.5,6.5,i);
lambdazVector= logspace(0.5,6.5,j);

path = [method,'_',eddy,'_',wall,'/singular_N=',num2str(N),'_R=',num2str(R),'_Nc=',num2str(Nc),'_i',num2str(i)];
load(path);

if strcmp(method,'SIOA')
    cormin = 0;
    if strcmp(eddy,'eddyoff')
        cormax = 4;
    else
        cormax = 2.5;
    end
else
    cormin = -4;
    if strcmp(eddy,'eddyoff')
        cormax = 6;
    else
        cormax = 2;
    end
end

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

A = log10(singular');

pcolor(lambdaxVector,lambdazVector,A);
set(gca,'xscale','log')
set(gca,'yscale','log')

xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');
shading interp

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
    if strcmp(eddy,'eddyoff')
    h.Ticks=[0,1,2,3,4];
    h.TickLabels={'0','1','2','3','4'};
    else
    h.Ticks=[0,0.5,1,1.5,2,2.5];
    h.TickLabels={'0','0.5','1','1.5','2','2.5'};
    end
end

set(gca,'xminortick','on'); 
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')

caxis([cormin,cormax]);
set(gcf,'unit','centimeters','position',[10 7 14 13.5]);
set(gca,'unit','centimeters','position',[3,3,9,9],'fontsize',fontn,'fontname','Times')


hold on
yline(100,'--k','LineWidth',2);
yline(0.6*R,'--k','LineWidth',2);
xline(1000,':k','LineWidth',2);
xline(6*R,':k','LineWidth',2);

% max point
[a,nxmax]=max(singular);
[c,nzmax]=max(a);
lambdaxmax = lambdaxVector(nxmax(nzmax));
lambdazmax = lambdazVector(nzmax);
fprintf('singularmax = %d \n',singular(nxmax(nzmax),nzmax));
fprintf('lambdaxmax = %d \n',lambdaxmax);
fprintf('lambdazmax = %d \n',lambdazmax);
loglog(lambdaxmax,lambdazmax,'^k','MarkerFaceColor','k','Markersize',10,'LineWidth',3)

% ridge line
if strcmp(method, 'SIOA')
    [sinmax,n] = max(singular); % fix lambda_z
    X = lambdaxVector(n(9:20));
    Z = lambdazVector(9:20); % lambda_z \in [30/R,3]
    K = polyfit(log10(X),log10(Z),1);
    k1 = K(1) % slope
    b1 = K(2);
    Xfit = (log10(Z)-b1)./k1;
    Xfit = 10.^Xfit;
end

print('-dpng','-r300', [path,'.png'])



