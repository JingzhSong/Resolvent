clc
clear
lambdax = 4000;
lambdaz = 890;
c = 16.6; 

R = 180;

kx = 2*pi/lambdax*R;
kz = 2*pi/lambdaz*R;

N = 100; % wall-normal grids

flowtype = 'turbu';
method = 'SIOA';
eddy = 'eddyoff';
wall = 'rigid'; % 'com', 'rigid'


i = sqrt(-1);
omega = c*kx;
T = 2*pi/omega;

phasemove =3;
fig1max=0.1; % amplitude limit
t0 = 0; % fix time for snap

Y=0;


[y,DM] = chebdif(N,2);
yplus = (y+1).*R;

fontn = 24;
ymax=180;

%% read
[u,v,w,p,f1,f2,f3,RSy,RS] = UF_up_w_new(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);

%% UVW & RS
wid = 36.5;
hei = 11.5;
figure;
set(gcf,'unit','centimeters','position',[2 10 wid hei])

subplot(1,4,1);
plot(abs(u),y,'k','LineWidth',2);hold on
plot(abs(v),y,'--r','LineWidth',2);hold on
plot(abs(w),y,':b','LineWidth',2);hold on
% plot(abs(p),y,'-m','LineWidth',1);hold on
grid on

axis([0,fig1max,-1,0]);
% legend('$\hat{u}$','$\hat{v}$','$\hat{w}$','$\hat{p}$','Interpreter','latex')
legend('$\hat{u}$','$\hat{v}$','$\hat{w}$','Interpreter','latex')
xlabel('$|\widehat{\mathbf u}|$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
% xticks(0:0.02:0.06)
set(gca,'unit','centimeters','position',[3,3,5,7],'fontsize',fontn,'fontname','Times');

subplot(1,4,2);
au = angle(u);
av = angle(v);
aw = angle(w);
ap = angle(p);
for j = 1:N
au(j) = au(j) + pi + phasemove;
if au(j) > 2 * pi
    au(j) = au(j) - 2*pi;
end
end
for j = 1:N
av(j) = av(j) + pi + phasemove;
if av(j) > 2 * pi
    av(j) = av(j) - 2*pi;
end
end
for j = 1:N
aw(j) = aw(j) + pi + phasemove;
if aw(j) > 2 * pi
    aw(j) = aw(j) - 2*pi;
end
end
for j = 1:N
ap(j) = ap(j) + pi + phasemove;
if ap(j) > 2 * pi
    ap(j) = ap(j) - 2*pi;
end
end
plot(au./2./pi,y,'k','LineWidth',2);hold on
plot(av./2./pi,y,'--r','LineWidth',2);hold on
plot(aw./2./pi,y,':b','LineWidth',2);hold on
% plot(ap./2./pi,y,'-m','LineWidth',1);hold on
grid on

axis([0,1,-1,0]);
% legend('$\hat{u}$','$\hat{v}$','$\hat{w}$','$\hat{p}$','location','NorthWest','Interpreter','latex')
xlabel('$(\angle{\widehat{\textbf{u}}})/2\pi$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'unit','centimeters','position',[11,3,5,7],'fontsize',fontn,'fontname','Times');

subplot(1,4,3);
plot(-1*real(RSy)*10^5,yplus,'k','LineWidth',2);hold on

grid on
% axis([-fig2max,fig2max,0,ymax]);
xlabel('$-{\rm Re}(\hat{u} ^ \ast \hat{v})$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'unit','centimeters','position',[19,3,5,7],'fontsize',fontn,'fontname','Times');
grid on

%%
subplot(1,4,4);
lambda_x = 2*pi/kx;
lambda_z = 2*pi/kz;
x = linspace(0,lambda_x,20);
z = linspace(0,lambda_z,20);

z0 = 0;
x0 = 0;
u_yz = real(u*exp(1i*(kx.*x0-omega.*t0)))*cos(kz.*z);
v_yz = real(v*exp(1i*(kx.*x0-omega.*t0)))*cos(kz.*z);
w_yz = -imag(w*exp(1i*(kx.*x0-omega.*t0)))*sin(kz.*z);

pcolor(z.*R,(y+1).*R,v_yz)
shading interp
axis([0,lambda_z*R,-inf,ymax]);
colormap(redblue);

h=colorbar;
h.Title.String = '${v}$';
h.Title.Interpreter = 'latex';
h.Title.FontSize = fontn;

set(gca,'Layer','top');
if strcmp(wall,'com')
    set(gca,'unit','centimeters','position',[27,2.6,6,7.45],'fontsize',fontn,'fontname','Times');
else
    set(gca,'unit','centimeters','position',[27,3,6,7],'fontsize',fontn,'fontname','Times');
end

xlabel('${z^+}$','Interpreter','latex');
ylabel('${y^+}$','Interpreter','latex');
grid on

hold on
quiver(z.*R,(y+1).*R,w_yz,v_yz,0.6,'Color','k','LineWidth',0.7);

% % plot wall displacement
% if strcmp(wall,'com')
% 
%     eta = 1i/omega*v(end,1);
% 
%     eta_yz = real(eta*exp(1i*(kx.*x0-omega.*t0)))*cos(kz.*z);
%     % eta_yz = real((eta*exp(1i*(kx.*x0+kz*z-omega.*t0))));  % wall displacement
%     hold on
%     plot(z.*R,eta_yz.*R*scaleeta,'k','LineWidth',2)
% end

%%
path = [method,'_',eddy,'_',wall,'/uvw_mushtaq2023_kx=',num2str(kx),'_kz=',num2str(kz),'_c=',num2str(c)];
print(gcf, [path,'.png'], '-dpng', '-r300');
