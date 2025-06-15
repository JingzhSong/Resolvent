clear
clc

N = 300;
R = 3300;
nc = 201;
nkx = 120;
nkz = 120;
flowtype = 'turbu';
method = 'IOA';
eddy = 'eddyoff';
wall = 'com';
Y = [];
uvwp = 'w';

path0 = [wall,'/',uvwp];
if strcmp(wall, 'rigid')
    path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),eddy,'kplus-70_final'];
    
else
    Cm = 0.46; Ck = 181; Cd=0.091;
    path1 = ['_N=',num2str(N),'nc=',num2str(nc),'nkx=',num2str(nkx),'nkz=',num2str(nkz),'_R=',num2str(R),'_Cm=',num2str(Cm),'_Ck=',num2str(Ck),'_Cd=',num2str(Cd),eddy,'kplus-70_final'];
end

% path = [path0,'/psd',path1];
% load([path,'.mat'],'psd');
path = [path0,'/cmax',path1];
load([path,'.mat'],'cmax');

kxVector = logspace(-7,0,nkx).*R;
kzVector = logspace(-7,0,nkz).*R;

lambdaxVector = 2.*pi./kxVector .*R;
lambdazVector = 2.*pi./kzVector .*R;

[y,~] = chebdif(N,2);
yplus = (y+1)*R;

if strcmp(wall, 'rigid')
    I = eye(N);
    k = 0.426;
    alpha = 25.4;
    NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
    nuT0 = NuT(y);
    DUDy = @(y) R.*(-y)./NuT(y);
    U1 = DUDy(y);
    U0 = zeros(N,1);
    for j=1:N
        U0(j) = integral(DUDy,-1,y(j));
    end
end
if strcmp(wall, 'com')
    if R==3300
        load('COMUmean3300.mat','Umean');
        U0 = Umean;
    end
    if R==6700
        load('COMUmean6700.mat','Umean');
        U0 = Umean;
    end
    if R==8900
        load('COMUmean8900.mat','Umean');
        U0 = Umean;
    end
end

%%
figure;
set(gcf,'unit','centimeters','position',[10 7 38 13]);

cormin = 1;
if strcmp(wall, 'com')
    cormax = 15;
else
    cormax =4;
end
fontn = 24;

subplot(1,3,1)
yi = length(find(yplus > 5));
fprintf('y=%d \n',yplus(yi));
c = cmax(:,:,yi);

% Divided by the local speed
pcolor(lambdaxVector,lambdazVector,c'./U0(yi)); hold on

set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+ $','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');

if R == 8900
    title('$ y^+ \approx 4 $','Interpreter','latex');
else
    title('$ y^+ \approx 7 $','Interpreter','latex');
end

shading interp
load('MyColormaps.mat')
colormap(mycmap);

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
axis([50 5*10^5 50 5*10^5])
caxis([cormin,cormax]);

% Draw large scale boundaries
hold on
xline(2*R,'--k','linewidth',2)
hold on
yline(0.4*R,'--k','linewidth',2)

set(gca,'unit','centimeters','position',[3.2,2.7,9,9],'fontsize',fontn,'fontname','Times')


%%
subplot(1,3,2)
yi = length(find(yplus > 15));
fprintf('y=%d \n',yplus(yi));
c = cmax(:,:,yi);

% Divided by the local speed
pcolor(lambdaxVector,lambdazVector,c'./U0(yi)); hold on

set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');


if R == 3300
    title('$ y^+ \approx 18 $','Interpreter','latex');
end
if R == 8900
    title('$ y^+ \approx 18 $','Interpreter','latex');
end
if R==6700
    title('$ y^+ \approx 17 $','Interpreter','latex');
end

shading interp
load('MyColormaps.mat')
colormap(mycmap);

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
axis([50 5*10^5 50 5*10^5])
caxis([cormin,cormax]);

% Draw large scale boundaries
hold on
xline(2*R,'--k','linewidth',2)
hold on
yline(0.4*R,'--k','linewidth',2)

set(gca,'unit','centimeters','position',[15,2.7,9,9],'fontsize',fontn,'fontname','Times')
set(gca,'fontsize',fontn,'fontname','Times')

%%
subplot(1,3,3)
yi = length(find(yplus > 105));
fprintf('y=%d \n',yplus(yi));
c = cmax(:,:,yi);

% Divided by the local speed
pcolor(lambdaxVector,lambdazVector,c'./U0(yi)); hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
xticks(logspace(1,6,6))
yticks(logspace(1,6,6))
xlabel('${\lambda}_x^+$','Interpreter','latex');
ylabel('${\lambda}_z^+$','Interpreter','latex');

if R == 3300
    title('$ y^+ \approx 113 $','Interpreter','latex');
else
    title('$ y^+ \approx 110 $','Interpreter','latex');
end

shading interp
load('MyColormaps.mat')
colormap(mycmap);
h=colorbar;
h=colorbar;
h.Label.Interpreter = 'latex';
h.Title.FontSize = fontn;
if strcmp(wall, 'com')
h.Ticks=[1,5,10,15,20];
else
h.Ticks=[1,2,3,4,5];
end

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.01]);
set(gca,'layer','top')
axis([50 5*10^5 50 5*10^5])
caxis([cormin,cormax]);

% Draw large scale boundaries
hold on
xline(2*R,'--k','linewidth',2)
hold on
yline(0.4*R,'--k','linewidth',2)

set(gca,'unit','centimeters','position',[26.7,2.7,9,9],'fontsize',fontn,'fontname','Times')

print('-dpng','-r300', [path0,'/cmaxk',path1,'.png'])





