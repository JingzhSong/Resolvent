clear
clc
R = 3300;

if R == 3300
    load('EXPUmean3300.mat')
    EXPUmean = EXPUmean3300;
    N = 300;
%     min_i = 6;
    min_i = 6;
    max_i = 76;
    uc = 26;
end
if R == 6700
    load('EXPUmean6700.mat')
    EXPUmean = EXPUmean6700;
    N = 400;
    min_i = 6;
    max_i = 75;
    uc = 28;
end
if R == 8900
    load('EXPUmean8900.mat')
    EXPUmean = EXPUmean8900;
    N = 400;
    min_i = 5;
    max_i = 65;
    uc = 28;
end
yplus0 = EXPUmean(:,1);
ypluslog0 = log10(EXPUmean(:,1));
Umean0 = EXPUmean(:,2);

[y,~] = chebdif(N,2);
yplus = (1-y(1:N/2,1)).*R;
ypluslog = log10(yplus);
yplusln = log(yplus);
Umean = zeros(N/2,1);
% Umean3300(6:76) = spline(ypluslog0,Umean0,ypluslog(6:76));

k=0.426;
alpha=25.4;
NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
DUDy = @(y) R.*(-y)./NuT(y);
U0 = zeros(N,1);
for j=1:N
    U0(j) = integral(DUDy,-1,y(j));
end

Umean = spline(ypluslog0,Umean0,ypluslog);
% semilogx(yplus,Umean(1:N/2),'linewidth',2)
% hold on
% semilogx(yplus,U0(1:N/2),'linewidth',2)
% hold on
% semilogx(yplus0,Umean0,'bo','linewidth',1)
% hold on
% grid on
% xlabel('y+')
% ylabel('U')
% legend('compliant','rigid','EXP&DNS')

% Umean = interp1(ypluslog0,Umean0,ypluslog);
% 内区u=y+c平移（对数坐标平移）
Umean(1:(min_i-1),1) = 10 .^ (ypluslog(1:(min_i-1),1) - ypluslog(min_i,1) + log10(Umean(min_i,1)));
% 内区u=y+c平移（非对数坐标平移）
% Umean(1:(min_i-1),1) = yplus(1:(min_i-1),1) + Umean(min_i,1) - yplus(min_i,1);
%外区：
% 对数区直接延长
% Umean((max_i+1):end,1) = 1./0.41 .*(yplusln(max_i+1):end,1) - yplusln(max_i,1)) + Umean(max_i,1);
% moser 文章（3.3）式
% Umean((max_i+1):end,1) = (1./0.402+150./R) .*(yplusln((max_i+1):end,1) - yplusln(max_i,1)) + 1/R .* (yplus((max_i+1):end,1) - yplus(max_i,1)) + Umean(max_i,1);
% moser文章（3.4）式
% Umean((max_i+1):end,1) = (1./0.384) .*(log(yplus((max_i+1):end,1)+1) - log(yplus(max_i,1)+1) ) + 1/2 * 0.7 / R^2 .* (yplus((max_i+1):end,1).^2 - yplus(max_i,1).^2) + Umean(max_i,1);
% 刚性壁涡黏模型直接平移
% Umean((max_i+1):end,1) = U0((max_i+1):N/2,1) + Umean(max_i,1) - U0(max_i,1);
% 刚性壁涡黏模型直接平移+峰值速度匹配
Umean((max_i+1):end,1) = (U0((max_i+1):N/2,1)- U0(max_i,1)) .* (uc - Umean(max_i,1)) ./ (U0(N/2,1) - U0(max_i,1)) + Umean(max_i,1)  ;

Umean = [Umean;flipud(Umean)];
Umean(Umean<0)=0;

U1 = gradient(U0,y);
U1mean = gradient(Umean,y);

% 用spline的解析式直接求导
% pp = spline(y,Umean);
% pp_deri = fnder(pp);
% U1mean = ppval(pp_deri,y);

semilogx(yplus,Umean(1:N/2),'linewidth',2)
hold on
semilogx(yplus,U0(1:N/2),'linewidth',2)
hold on
semilogx(yplus0,Umean0,'bo','linewidth',1)
hold on
grid on
xlabel('y+')
ylabel('U')
legend('compliant','rigid','EXP&DNS')

figure;
plot(y,Umean,'linewidth',1)
hold on
plot(y,U0,'linewidth',1)
grid on
xlabel('y')
ylabel('U')
legend('compliant','rigid')

% figure;
% semilogx(yplus(N/2:end,1),Umean(N/2:end,1),'linewidth',2);hold on
% semilogx(yplus(N/2:end,1),U0(N/2:end,1),'linewidth',2);
% xlabel('y+')
% ylabel('U')
% legend('rigid','compliant')

figure;
yplus = (y+1)*R;
semilogx(yplus(N/2:end,1),U1mean(N/2:end,1),'linewidth',2); hold on
semilogx(yplus(N/2:end,1),U1(N/2:end,1),'linewidth',2); 
xlabel('y+')
ylabel('dU/dy')
legend('compliant','rigid')

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
N = 300;
R = 3300;

min_i = 6;
max_i = 75;
EXPUmean = EXPUmean3300;

yplus0 = EXPUmean(:,1);

ypluslog0 = log10(EXPUmean(:,1));

Umean0 = EXPUmean(:,2);

[y,~] = chebdif(N,2);
yplus = (1-y(1:N/2,1)).*R;
ypluslog = log10((1-y(1:N/2,1)).*R);
yplusln = log((1-y(1:N/2,1)).*R);
Umean = zeros(N/2,1);
% Umean3300(6:76) = spline(ypluslog0,Umean0,ypluslog(6:76));

Umean = interp1(ypluslog0,Umean0,ypluslog);
Umean(1:(min_i-1),1) = 10 .^ (ypluslog(1:(min_i-1),1) - ypluslog(min_i,1) + log10(Umean(min_i,1)));
Umean(76:end,1) = 1./0.41 .*(yplusln(76:end,1) - yplusln(75,1)) + Umean(75,1);
% Umean((max_i+1):end,1) = (1./0.402+150./2000) .*(yplusln((max_i+1):end,1) - yplusln(max_i,1)) + 1/2000 .* (yplus((max_i+1):end,1) - yplus(max_i,1)) + Umean(max_i,1);

Umean = [Umean;flipud(Umean)];

% plot(ypluslog,Umean6700);hold on
% plot(ypluslog0,Umean0);hold on

% [y,~] = chebdif(N,2);
%     yplus = (y+1)*R;
%     I = eye(N);
%     k=0.426;
%     alpha=25.4;
%     NuT = @(y) 0.5.*(1+(k.*R./3.*(2.*(y+1)-(y+1).^2).*(3-4.*(y+1)+2.*(y+1).^2).*(1-exp((abs(y)-1).*R./alpha))).^2).^0.5 + 0.5;
%     DUDy = @(y) R.*(-y)./NuT(y);
%     U0 = zeros(N,1);
%     for j=1:N
%         U0(j) = integral(DUDy,-1,y(j));
%     end
%
% plot(y,U0)

semilogx(yplus,Umean(1:N/2))
hold on
semilogx(yplus0,Umean0,'bo')
hold on

%% 取点正确
ypluslog0 = log10(EXPUmean3300(:,1));
semilogx(EXPUmean3300(:,1),EXPUmean3300(:,2));hold on
ypluslog0 = log10(EXPUmean6700(:,1));
semilogx(EXPUmean6700(:,1),EXPUmean6700(:,2));hold on
ypluslog0 = log10(EXPUmean8900(:,1));
semilogx(EXPUmean8900(:,1),EXPUmean8900(:,2));hold on
legend('3300','6700','8900')

%%

% semilogx(EXPUmean3300(:,1)./64,EXPUmean3300(:,2));hold on
% Uc = trapz(EXPUmean3300(:,1),EXPUmean3300(:,2)) ./ (EXPUmean3300(end,1) - EXPUmean3300(1,1));

load('COMUmean3300.mat','Umean');
U0 = Umean;
% Uc = U0(N/2);
N=300;R=3300;
[y,~] = chebdif(N,2);
yplus = (y+1)*R;

% Uc = -trapz(y(1:N/2),Umean3300(1:N/2)) ;
Uc = 26;
semilogx(yplus(N/2:end,1)./64,U0(N/2:end,1)./Uc,'linewidth',2);hold on
semilogx(EXPu3300(:,1),EXPu3300(:,2),'o','linewidth',2);hold on


%%

% semilogx(EXPUmean6700(:,1)./165,EXPUmean6700(:,2));hold on
% Uc = trapz(EXPUmean6700(:,1),EXPUmean6700(:,2)) ./ (EXPUmean6700(end,1) - EXPUmean6700(1,1));

load('COMUmean6700.mat','Umean');
U0 = Umean;
% Uc = U0(N/2);

N=400;R=6700;
[y,~] = chebdif(N,2);
yplus = (y+1)*R;
% Uc = -trapz(y(1:N/2),Umean6700(1:N/2)) ;

Uc=28;
semilogx(yplus(N/2:end,1)./165,U0(N/2:end,1)./Uc,'linewidth',2);hold on
% semilogx(EXPu6700(:,1),EXPu6700(:,2),'o','linewidth',2);hold on


%%

% semilogx(EXPUmean8900(:,1)./192,EXPUmean8900(:,2));hold on
% Uc = trapz(EXPUmean8900(:,1),EXPUmean8900(:,2)) ./ (EXPUmean8900(end,1) - EXPUmean8900(1,1));

load('COMUmean8900.mat','Umean');
U0 = Umean;
% Uc = U0(N/2);
N=400;R=8900;
[y,~] = chebdif(N,2);
yplus = (y+1)*R;
% Uc = -trapz(y(1:N/2),Umean8900(1:N/2)) ;

Uc = 28;
semilogx(yplus(N/2:end,1)./192,U0(N/2:end,1)./Uc,'linewidth',2);hold on
semilogx(EXPu8900(:,1),EXPu8900(:,2),'o','linewidth',2);hold on

legend('Re3300 Umean','Re3300 Ucov','Re6700 Umean','Re6700 Ucov','Re8900 Umean','Re8900 Ucov')
xlabel('y+/yc')
ylabel('U/U0')

