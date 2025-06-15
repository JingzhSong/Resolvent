clear
clc
kx = 12;kz = 120;c = 10;  % near-wall cycle
% kx = 1;kz = 10;c = 16;   % VLSMs
omega = c*kx;

N = 200;  % wall-normal grid
NY = 500; % Cd and Ck grid
R = 2000;  %Renolds number
i = sqrt(-1);

cormin=-1;
cormax=1;

fontn = 24;

Cm = 2;
if kx == 12
    Ckvector = linspace(27000,31000,NY);
    Cdvector = linspace(-10,10,NY);
end
if kx ==1
    Ckvector = linspace(350,700,NY);
    Cdvector = linspace(-10,10,NY);
end
[Cd,Ck] = meshgrid(Cdvector,Ckvector);
[ReY,ImY] = C2Y(Cm,Ck,Cd,omega);

pcolor(Cd,Ck,ReY)
colormap(redblue);
caxis([cormin,cormax]);
colorbar
shading interp
hold on
Ckfun = @(Cdvar) (2.*Cm.^2.*omega.^2 + Cdvar.^2) ./ 2 ./ Cm;
Ckr = Ckfun(Cdvector);
plot(Cdvector,Ckr,'k','LineWidth',1)
xlabel('${C_d}$','Interpreter','latex');
ylabel('${C_k}$','Interpreter','latex');
h=colorbar;
h.Title.String = '${\rm Re}(Y)$';
h.Title.Interpreter = 'latex';
h.Title.FontSize = fontn;

 set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.015]);
set(gca,'layer','top');
set(gcf,'unit','centimeters','position',[10 7 15 12.5]);
set(gca,'unit','centimeters','position',[3,2.2,9,9],'fontsize',fontn,'fontname','Times')

print('-dpng','-r300', ['ReY_kx12.png'])

figure;
pcolor(Cd,Ck,ImY)
colormap(redblue);
caxis([cormin,cormax]);
colorbar
shading interp
hold on

Ckfun = @(Cdvar) (2.*Cm.^2.*omega.^2 + Cdvar.^2) ./ 2 ./ Cm;
Ckr = Ckfun(Cdvector);
plot(Cdvector,Ckr,'k','LineWidth',1)

xlabel('${C_d}$','Interpreter','latex');
ylabel('${C_k}$','Interpreter','latex');
h=colorbar;
h.Title.String = '${\rm Im}(Y)$';
h.Title.Interpreter = 'latex';
h.Title.FontSize = fontn;

set(gca,'xminortick','on');
set(gca,'ticklength',[0.02 0.015]);
set(gca,'layer','top');

set(gcf,'unit','centimeters','position',[10 7 15 12.5]);
set(gca,'unit','centimeters','position',[3,2.2,9,9],'fontsize',fontn,'fontname','Times')

print('-dpng','-r300', ['ImY_kx12.png'])