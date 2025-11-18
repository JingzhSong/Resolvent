%% 基本参数等
clc
clear
% kx = 12;kz = 120;c = 10;  %Fig.5(a)
% kx = 1;kz = 10;c = 16;  %Fig.7(a)
kx = 0.0675;kz = 0.1422;c = 0;
N = 60;
R = 350;  %Renolds number

flowtype = 'couette';
method = 'SIOA';
eddy = 'eddyoff';
wall = 'rigid'; % 'com', 'rigid'

omega = c*kx;
i = sqrt(-1);

[y,DM] = chebdif(N,2);
D1 = DM(1:N,1:N,1);
D2 = DM(1:N,1:N,2);
I = eye(N);
grd = [i*kx.*I; D1; i*kz.*I];
lplc = D2 - kx.^2.*I - kz.^2.*I;

%% read
Y=0;

[u,v,w,p,f1,f2,f3,RSy,RS] = UF_up_w_new(kx,kz,c,N,R,flowtype,method,eddy,wall,Y);

%%
fontn = 22;

[y,DM] = chebdif(N,2);

figure;
set(gcf,'unit','centimeters','position',[2 10 35 12])

subplot(1,3,1);
plot(y,real(u),'--r','LineWidth',2);hold on
plot(y,imag(u),':b','LineWidth',2);hold on
plot(y,abs(u),'-k','LineWidth',2);hold on
axis([-1,1,-0.1,0.4]);
legend('${\rm Re}(u)$','${\rm Im}(u)$','$|u|$','Interpreter','latex')
% xlabel('$\rm Re(u)$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'fontsize',fontn,'fontname','Times');
grid on;

subplot(1,3,2);
plot(y,real(v),'--r','LineWidth',2);hold on
plot(y,imag(v),':b','LineWidth',2);hold on
plot(y,abs(v),'-k','LineWidth',2);hold on
axis([-1,1,-0.02,0.02]);
legend('${\rm Re}(v)$','${\rm Im}(v)$','$|v|$','Interpreter','latex')
% xlabel('$\rm Re(v)$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'fontsize',fontn,'fontname','Times');
grid on;

subplot(1,3,3);
plot(y,real(w),'--r','LineWidth',2);hold on
plot(y,imag(w),':b','LineWidth',2);hold on
plot(y,abs(w),'-k','LineWidth',2);hold on
axis([-1,1,-0.3,0.45]);
legend('${\rm Re}(w)$','${\rm Im}(w)$','$|w|$','Interpreter','latex')
% xlabel('$\rm Re(w)$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'fontsize',fontn,'fontname','Times');
grid on;


%%
fontn = 22;

[y,DM] = chebdif(N,2);

u=f1;
v=f2;
w=f3;
figure;
set(gcf,'unit','centimeters','position',[2 10 35 12])

subplot(1,3,1);
plot(y,real(u),'--r','LineWidth',2);hold on
plot(y,imag(u),':b','LineWidth',2);hold on
plot(y,abs(u),'-k','LineWidth',2);hold on
% axis([-1,1,-0.1,0.4]);
legend('${\rm Re}(f_x)$','${\rm Im}(f_x)$','$|f_x|$','Interpreter','latex')
% xlabel('$\rm Re(u)$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'fontsize',fontn,'fontname','Times');
grid on;

subplot(1,3,2);
plot(y,real(v),'--r','LineWidth',2);hold on
plot(y,imag(v),':b','LineWidth',2);hold on
plot(y,abs(v),'-k','LineWidth',2);hold on
% axis([-1,1,-0.02,0.02]);
legend('${\rm Re}(f_y)$','${\rm Im}(f_y)$','$|f_y|$','Interpreter','latex')
% xlabel('$\rm Re(v)$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'fontsize',fontn,'fontname','Times');
grid on;

subplot(1,3,3);
plot(y,real(w),'--r','LineWidth',2);hold on
plot(y,imag(w),':b','LineWidth',2);hold on
plot(y,abs(w),'-k','LineWidth',2);hold on
% axis([-1,1,-0.3,0.45]);
legend('${\rm Re}(f_z)$','${\rm Im}(f_z)$','$|f_z|$','Interpreter','latex')
% xlabel('$\rm Re(w)$','Interpreter','latex');
ylabel('$y^+$','Interpreter','latex');
set(gca,'fontsize',fontn,'fontname','Times');
grid on;