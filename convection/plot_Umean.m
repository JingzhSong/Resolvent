clear
clc

R = 3300;
load('EXPUmean3300.mat')
EXPUmean = EXPUmean3300;
N = 300;
min_i = 6;
max_i = 76;
uc = 26;

yplus0 = EXPUmean(:,1);
ypluslog0 = log10(EXPUmean(:,1));
Umean0 = EXPUmean(:,2);

[y,~] = chebdif(N,2);
yplus = (1-y(1:N/2,1)).*R;
ypluslog = log10(yplus);
yplusln = log(yplus);
Umean = zeros(N/2,1);

k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
DUDy = @(y) R.*(-y)./NuT(y);
U0 = zeros(N,1);
for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end

Umean = interp1(ypluslog0,Umean0,ypluslog(1:N/2));
% inner region: u=y+c in log axis
Umean(1:(min_i-1),1) = 10 .^ (ypluslog(1:(min_i-1),1) - ypluslog(min_i,1) + log10(Umean(min_i,1)));
% outer region: rigid wall eddy model - Translation + peak velocity matching scaling
Umean((max_i+1):end,1) = (U0((max_i+1):N/2,1)- U0(max_i,1)) .* (uc - Umean(max_i,1)) ./ (U0(N/2,1) - U0(max_i,1)) + Umean(max_i,1)  ;

Umean = [Umean;flipud(Umean)];
Umean(Umean<0)=0;

U0_3300 = U0;
Umean0_3300 = Umean0;
yplus_3300 = yplus;
yplus0_3300 = yplus0;

if R ==3300
    save('COMUmean3300.mat','Umean')
end
if R ==6700
    save('COMUmean6700.mat','Umean')
end
if R ==8900
    save('COMUmean8900.mat','Umean')
end

%%
R = 6700;

load('EXPUmean6700.mat')
EXPUmean = EXPUmean6700;
N = 400;
min_i = 6;
max_i = 75;
uc = 28;

yplus0 = EXPUmean(:,1);
ypluslog0 = log10(EXPUmean(:,1));
Umean0 = EXPUmean(:,2);

[y,~] = chebdif(N,2);
yplus = (1-y(1:N/2,1)).*R;
ypluslog = log10(yplus);
yplusln = log(yplus);
Umean = zeros(N/2,1);

k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
DUDy = @(y) R.*(-y)./NuT(y);
U0 = zeros(N,1);
for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end

Umean = interp1(ypluslog0,Umean0,ypluslog(1:N/2));
% inner region: u=y+c in log axis
Umean(1:(min_i-1),1) = 10 .^ (ypluslog(1:(min_i-1),1) - ypluslog(min_i,1) + log10(Umean(min_i,1)));
% outer region: rigid wall eddy model - Translation + peak velocity matching scaling
Umean((max_i+1):end,1) = (U0((max_i+1):N/2,1)- U0(max_i,1)) .* (uc - Umean(max_i,1)) ./ (U0(N/2,1) - U0(max_i,1)) + Umean(max_i,1)  ;

Umean = [Umean;flipud(Umean)];
Umean(Umean<0)=0;

U0_6700 = U0;
Umean0_6700 = Umean0;
yplus_6700 = yplus;
yplus0_6700 = yplus0;

if R ==3300
    save('COMUmean3300.mat','Umean')
end
if R ==6700
    save('COMUmean6700.mat','Umean')
end
if R ==8900
    save('COMUmean8900.mat','Umean')
end

%%
R = 8900;

load('EXPUmean8900.mat')
EXPUmean = EXPUmean8900;
N = 400;
min_i = 5;
max_i = 65;
uc = 28;

yplus0 = EXPUmean(:,1);
ypluslog0 = log10(EXPUmean(:,1));
Umean0 = EXPUmean(:,2);

[y,~] = chebdif(N,2);
yplus = (1-y(1:N/2,1)).*R;
ypluslog = log10(yplus);
yplusln = log(yplus);
Umean = zeros(N/2,1);

k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
DUDy = @(y) R.*(-y)./NuT(y);
U0 = zeros(N,1);
for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end

Umean = interp1(ypluslog0,Umean0,ypluslog(1:N/2));
% inner region: u=y+c in log axis
Umean(1:(min_i-1),1) = 10 .^ (ypluslog(1:(min_i-1),1) - ypluslog(min_i,1) + log10(Umean(min_i,1)));
% outer region: rigid wall eddy model - Translation + peak velocity matching scaling
Umean((max_i+1):end,1) = (U0((max_i+1):N/2,1)- U0(max_i,1)) .* (uc - Umean(max_i,1)) ./ (U0(N/2,1) - U0(max_i,1)) + Umean(max_i,1)  ;

Umean = [Umean;flipud(Umean)];
Umean(Umean<0)=0;

Umean_8900 = Umean;
U0_8900 = U0;
Umean0_8900 = Umean0;
yplus_8900 = yplus;
yplus0_8900 = yplus0;

if R ==3300
    save('COMUmean3300.mat','Umean')
end
if R ==6700
    save('COMUmean6700.mat','Umean')
end
if R ==8900
    save('COMUmean8900.mat','Umean')
end

%%
semilogx(0,0,'w')
hold on
load('COMUmean3300.mat','Umean')
Umean_3300 = Umean;
load('COMUmean6700.mat','Umean')
Umean_6700 = Umean;
load('COMUmean8900.mat','Umean')
Umean_8900 = Umean;
semilogx(yplus_3300,Umean_3300(1:150),'k-','linewidth',2)
hold on
semilogx(yplus_6700,Umean_6700(1:200),'b-','linewidth',2)
hold on
semilogx(yplus_8900,Umean_8900(1:200),'r-','linewidth',2)
hold on

semilogx(0,0,'w')
hold on
semilogx(yplus_3300,U0_3300(1:150),'k--','linewidth',2)
hold on
semilogx(yplus_6700,U0_6700(1:200),'b--','linewidth',2)
hold on
semilogx(yplus_8900,U0_8900(1:200),'r--','linewidth',2)
hold on

nsize = 12;
semilogx(0,0,'w')
hold on
semilogx(yplus0_3300(7:end),Umean0_3300(7:end),'k+','linewidth',2,'Markersize',nsize)
hold on
semilogx(yplus0_6700(7:end),Umean0_6700(7:end),'bo','linewidth',2,'Markersize',nsize)
hold on
semilogx(yplus0_8900(8:end),Umean0_8900(8:end),'r*','linewidth',2,'Markersize',nsize)
hold on

semilogx(0,0,'w')
hold on
h = semilogx(yplus0_3300(1:6),Umean0_3300(1:6),'m^','linewidth',2,'Markersize',nsize);
hold on

grid on
fontn = 22;
xlabel('$y^+$','Interpreter','latex');
ylabel('$U^+$','Interpreter','latex');
set(gcf,'unit','centimeters','position',[10 7 20 13]);
set(gca,'unit','centimeters','position',[3.5,2.5,15,9],'fontsize',fontn,'fontname','Times')

print('-dpng','-r300', ['Umean.png'])

legend('Compliant wall','$Re_{\tau}=3300$','$Re_{\tau}=6700$','$Re_{\tau}=8900$','Rigid wall','$Re_{\tau}=3300$','$Re_{\tau}=6700$','$Re_{\tau}=8900$','Lu $et~al.~(2024)$','$Re_{\tau}=3300$','$Re_{\tau}=6700$','$Re_{\tau}=8900$','Wang $et~al.~(2020)$','$Re_{\tau}=8600$','Interpreter','latex','location','eastoutside','fontsize',14,'NumColumns', 1)
print('-dpng','-r300', ['Umean_colorbar.png'])
