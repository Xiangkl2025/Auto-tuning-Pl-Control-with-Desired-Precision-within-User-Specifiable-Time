clc; 
clear all ;
step=0.0001;
n=200000;
t=0:step:step*n;
%% design parameters
gamma=0.5;
K_D=2;
n1=3.5;n2=0.76;n3=0.8;n4=0.3;n5=0.5;n6=0.5;
T=3;k1=5;k2=5;
%% initial value
x1=0.5;
x2=0.6;
dx1=0.3;
dx2=0.6;
ddx2=0.3;
W1=1;W2=1;
delta1=0.1;
delta2=0.1;
rhof1=0.1;
rhof2=0.1;
beta1=0.6;
beta2=0.4;
c=1;
a1=0.1;a2=20;Y1=0;Y2=0;b1=30;b2=20;sigma=0.5;
hat_a=0;
%%
for i=1:1:n;
    tt=i*step;
%%
x1s(i)=x1;
x2s(i)=x2;
%% desired signals and its derivative
yd1=0.3*cos(0.5*tt); 
dyd1=-0.3*0.5*sin(0.5*tt);
ddyd1=-0.3*0.5*0.5*cos(0.5*tt);
yd1s(i)=yd1;
dyd1s(i)=dyd1;
ddyd1s(i)=ddyd1;
yd2=0.4*cos(0.5*tt); 
dyd2=-0.4*0.5*sin(0.5*tt);
ddyd2=-0.4*0.5*0.5*cos(0.5*tt);
yd2s(i)=yd2;
dyd2s(i)=dyd2;
ddyd2s(i)=ddyd2;
%% tracking error 
z1=x1-yd1;
z1s(i)=z1;
z2=x2-yd2;
z2s(i)=z2;
%% dz1
dz1=dx1-dyd1;
dz1s(i)=dz1;
ddz1=dx2-ddyd1;
ddz1s(i)=ddz1;
dz1=dx1-dyd1;
dz1s(i)=dz1;
ddz1=dx2-ddyd1;
ddz1s(i)=ddz1;
%% dz2
dz2=dx2-dyd2;
dz2s(i)=dz2;
ddz2=dx2-ddyd2;
ddz2s(i)=ddz2;
dz2=dx2-dyd2;
dz2s(i)=dz2;
ddz2=dx2-ddyd2;
ddz2s(i)=ddz2;
%% filiter error
Z1=W1*z1+dz1;
Z2=W2*z2+dz2;
%% rho
if tt < T  
   rho1=(1/tanh(k1*tt))*(T-tt)^2+rhof1;  
else  
   rho1=rhof1;  
    end 
if tt < T  
   rho2=(1/tanh(k2*tt))*(T-tt)^2+rhof2;  
else  
   rho2=rhof2;  
end 
   rho1s(i)=rho1;
    rho2s(i)=rho2;
    
      %% drho
if tt < T  
   drho1=-(1/(sinh(k1*tt)))^2*(-k1)*(T-tt)^2-3*(1/tanh(k1*tt))*(T-tt)^2;
else  
   drho1=0;  
    end 
if tt < T  
   drho2=-(1/(sinh(k2*tt)))^2*(-k2)*(T-tt)^2-3*(1/tanh(k2*tt))*(T-tt)^2;
else  
   drho2=0;  
end 
%% transformation error
X1=delta1*(rho1+Z1)*(rho1-Z1);
X2=delta2*(rho2+Z2)*(rho2-Z2);
%% h
if 0<X1 <c 
  h1=1-((X1/c)-1)^3;  
else  
  h1=1;  
    end 
if 0<X2 <c 
  h2=1-((X2/c)-1)^3;  
else  
  h2=1;  
end
%% varepsilon
varepsilon1=Z1/h1;
varepsilon2=Z2/h2;
%% r
if 0<X1 <c 
  r1=1/h1-(12*delta1/(c*h1^2))*(X1/c-1)^2*z1^2;  
else  
  r1=1;  
end 
 if 0<X2 <c 
  r2=1/h2-(12*delta2/(c*h2^2))*(X2/c-1)^2*z2^2;  
else  
  r2=1;  
 end  
    %% v
    if 0<X1 <c 
  v1=(12*delta1/(c*h1^2))*(X1/c-1)^2*rho1*drho1;  
else  
  v1=0;  
end 
 if 0<X2 <c 
  v2=(12*delta2/(c*h2^2))*(X2/c-1)^2*rho2*drho2;  
else  
  v2=0;  
 end   
 %%
R11=-rho1;R11s(i)=R11;
R12=rho1;R12s(i)=R12;
R21=-(rho2);R21s(i)=R21;
R22=rho2;R22s(i)=R22;
%% 积分 varepsilon1
Y1=Y1+(varepsilon1)*step;
Y2=Y2+(varepsilon2)*step;
%% 积分
Q1=Y1/(b1*Y1^1+1);
Q2=Y2/(b2*Y2^1+1);
Q1s(i)=Q1;
Q2s(i)=Q2;
%% dotH, dotG dotF
g1=beta1*(1/((a1*varepsilon1)^2+1))+beta2;
g2=beta1*(1/((a2*varepsilon2)^2+1))+beta2;
h1=(1-b1*Y1^2)/((b1*Y1)+1)^2;
h2=(1-b2*Y2^2)/((b2*Y2)+1)^2;
dotG=[g1;g2];dotH=[h1;h2];
g1s(i)=g1;g2s(i)=g2;
%% phi
varphi=3+x2^2;
varphis(i)=varphi;
%% E_1
  E_1=beta1*(varepsilon1/(sqrt((a1*varepsilon1)^2+1)))+beta2*varepsilon1+Y1/((b1*Y1)^2+1);
   E_1s(i)=E_1;
   %% E2
  E_2=beta2*(varepsilon2/(sqrt((a2*varepsilon2)^2+1)))+beta2*varepsilon2+Y2/((b2*Y2)^2+1);
   E_2s(i)=E_2;
   E=[E_1;E_2];
%% psi
psi=norm(dotG)*varphi+norm(dotH);
psis(i)=psi;
%% adaptive law   
hat_a=hat_a+step*(gamma*(norm(E)^2*psi^2-sigma*hat_a));
hat_as(i)=hat_a;
%% Control input u
DeltaK_D=hat_a*psi^2;
DeltaK_Ds(i)=DeltaK_D;
%% control input u
u1=-r1*g1*(K_D*E_1+DeltaK_D*E_1);
u2=-r2*g2*(K_D*E_2+DeltaK_D*E_2);
u1s(i)=u1;u2s(i)=u2;
%% plant
x1=x1+step*(dx1);
x2=x2+step*(dx2);
x1s(i)=x1; 
x2s(i)=x2; 
% H
h11=n1+n2+2*n3*cos(x2);
h12=n2+n3*cos(x2);
h21=n2+n3*cos(x2);
h22=n2;
% Ng
ng11=-n3*dx2*sin(x1);
ng12=-n3*(dx1+dx2)*sin(x2);
ng21=n3*dx1*sin(x2);
ng22=0;
% G
gg1=n4*9.8*cos(x1)+n5*9.8*cos(x1+x2)+sin(n6*tt);
gg2=n5*9.8*cos(x1+x2)+cos(n6*tt);
%d
d11=sin(0.5*tt)+5*dx1^7;
d12=cos(0.5*tt)+5*dx2^7;
ddx1=(-h12*ddx2-ng11*dx1-ng12*dx2-gg1-d11+u1)/h11;
ddx2=(-h21*ddx1-ng21*dx1-ng22*dx2-gg2-d12+u2)/h22;
ddx1s(i)=ddx1;
ddx2s(i)=ddx2;
% dq
dx1=dx1+step*(ddx1);
dx2=dx2+step*(ddx2);
dx1s(i)=dx1;  
dx2s(i)=dx2; 
end
%% 对比线性PI
gamma1=0.5;c=1;
K_D1=2;
n11=3.5;n21=0.76;n31=0.8;n41=0.3;n51=0.5;n61=0.5;
T1=4;k11=5;k21=5;
%% initial value
x11=0.5;
hat_a1=0;
x21=0.6;
dx11=0.3;
dx21=0.6;
ddx21=0.3;
W11=1;W21=1;
delta11=0.1;
delta21=0.1;
rhof11=0.1;
rhof21=0.1;
beta11=0.6;
beta21=0.4;
c1=1;
Y11=0;Y21=0;
sigma1=0.5;
%%
for i=1:1:n;
    tt=i*step;
%%
x11s(i)=x11;
x21s(i)=x21;
%% desired signals and its derivative
yd11=0.3*cos(0.5*tt); 
dyd11=-0.3*0.5*sin(0.5*tt);
ddyd11=-0.3*0.5*0.5*cos(0.5*tt);
yd11s(i)=yd11;
dyd11s(i)=dyd11;
ddyd11s(i)=ddyd11;
yd21=0.4*cos(0.5*tt); 
dyd21=-0.4*0.5*sin(0.5*tt);
ddyd21=-0.4*0.5*0.5*cos(0.5*tt);
yd21s(i)=yd21;
dyd21s(i)=dyd21;
ddyd21s(i)=ddyd21;
%% tracking error 
z111=x11-yd11;
z111s(i)=z111;
z222=x21-yd21;
z222s(i)=z222;z=[z111,z222];
%% dz1
dz111=dx11-dyd11;
dz111s(i)=dz111;
ddz111=dx21-ddyd11;
ddz111s(i)=ddz111;

%% dz2
dz222=dx21-dyd21;
dz222s(i)=dz222;
ddz222=dx21-ddyd21;
ddz222s(i)=ddz222;
%% phi
varphi11=3.5+x21^2;
varphi11s(i)=varphi11;
%% filiter error
Z11=W11*z111+dz111;
Z21=W21*z222+dz222;
%% rho
if tt < T1  
   rho11=(1/tanh(k11*tt))*(T1-tt)^2+rhof11;  
else  
   rho11=rhof11;  
    end 
if tt < T1  
   rho21=(1/tanh(k21*tt))*(T1-tt)^2+rhof21;  
else  
   rho21=rhof21;  
end 
   rho11s(i)=rho11;
    rho21s(i)=rho21;
    
      %% drho
if tt < T1  
   drho11=-(1/(sinh(k11*tt)))^2*(-k11)*(T1-tt)^2-3*(1/tanh(k11*tt))*(T1-tt)^2;
else  
   drho11=0;  
    end 
if tt < T1  
   drho21=-(1/(sinh(k21*tt)))^2*(-k21)*(T1-tt)^2-3*(1/tanh(k21*tt))*(T1-tt)^2;
else  
   drho21=0;  
end 
%% transformation error
X11=delta11*(rho11+Z11)*(rho11-Z11);
X21=delta21*(rho21+Z21)*(rho21-Z21);
%% h
if 0<X11 <c1
  h11=1-((X11/c1)-1)^3;  
else  
  h11=1;  
    end 
if 0<X21 <c 
  h21=1-((X21/c1)-1)^3;  
else  
  h21=1;  
end
%% varepsilon
varepsilon11=Z11/h11;
varepsilon21=Z21/h21;
%% r
if 0<X11 <c 
  r11=1/h11-(12*delta11/(c1*h11^2))*(X11/c1-1)^2*z111^2;  
else  
  r11=1;  
end 
 if 0<X21 <c 
  r21=1/h21-(12*delta21/(c1*h21^2))*(X21/c1-1)^2*z222^2;  
else  
  r21=1;  
 end  
    %% v
    if 0<X11 <c 
  v11=(12*delta11/(c*h11^2))*(X11/c1-1)^2*rho11*drho11;  
else  
  v11=0;  
end 
 if 0<X21 <c 
  v21=(12*delta21/(c1*h21^2))*(X21/c1-1)^2*rho21*drho21;  
else  
  v21=0;  
 end   
%% 积分 varepsilon1
Y11=Y11+(varepsilon11)*step;
Y21=Y21+(varepsilon21)*step;
Y11s(i)=Y11;
Y21s(i)=Y21;
%% E_11
  E_111=varepsilon11+Y11;
   E_111s(i)=E_111;
   %% E2
E_222=varepsilon21+Y21;
   E_222s(i)=E_222; 
   E11=[E_111;E_222];
   dz=[dz111;dz222];ddz=[ddz111;ddz222];

%% psi
psi11=varphi11+norm(Z11);
psi11s(i)=psi11;
%% adaptive law   
hat_a1=hat_a1+step*(gamma1*(norm(E11)^2*psi11^2)-hat_a1*sigma1);
hat_a1s(i)=hat_a1;
%% Control input u
DeltaK_D1=hat_a1*psi11^2;
DeltaK_D1s(i)=DeltaK_D1;
%% control input u
u111=-r11*(K_D1*E_111+DeltaK_D1*E_111);
u222=-r21*(K_D1*E_222+DeltaK_D1*E_222);
u111s(i)=u111;u222s(i)=u222; U11=(u111^2+u222^2)*step;

%% plant
x11=x11+step*(dx11);
x21=x21+step*(dx21);
x11s(i)=x11; 
x21s(i)=x21; 
% H
h111=n11+n21+2*n31*cos(x21);
h121=n21+n31*cos(x21);
h211=n21+n31*cos(x21);
h221=n21;
% Ng
ng111=-n31*dx21*sin(x21);
ng121=-n31*(dx11+dx21)*sin(x21);
ng211=n31*dx11*sin(x21);
ng221=0;
% G
gg11=n41*9.8*cos(x11)+n51*9.8*cos(x11+x21)+sin(n61*tt);
gg21=n51*9.8*cos(x11+x21)+cos(n61*tt);
%d
d111=sin(0.5*tt)+5*dx11^7;
d121=cos(0.5*tt)+5*dx21^7;
ddx11=(-h121*ddx21-ng111*dx11-ng121*dx21-gg11-d111+u111)/h111;
ddx21=(-h211*ddx11-ng211*dx11-ng221*dx21-gg21-d121+u222)/h221;
ddx11s(i)=ddx11;
ddx21s(i)=ddx21;
% dq
dx11=dx11+step*(ddx11);
dx21=dx21+step*(ddx21);
dx11s(i)=dx11;  
dx21s(i)=dx21; 
end
%%
 figure(1);
plot(t(1:n),z1s,'r',t(1:n),z111s,'b-',t(1:n),R11s,'m--',t(1:n),R12s,'m--','linewidth',2.5);xlabel('Time(sec)','FontName','Times New Roman');
legend('$e_1$ (The proposed nonlinear adaptive PI)','$e_1$ (Linear PI with time-varying gains by [15])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('$e_1$(rad) ','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
  ylim([-10 10]);
  hold on
plot([T T], ylim, 'k--', 'LineWidth', 2);
legend('$e_1$ (The proposed nonlinear adaptive PI)','$e_1$ (Linear PI with time-varying gains by [15])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
axes('position',[0.55 0.15 0.35 0.25]);
plot(t(1:n),z1s,'r',t(1:n),z111s,'b-',t(1:n),R11s,'m--',t(1:n),R12s,'m--','linewidth',2.5);
xlim([2 20]);
ylim([-0.1 0.1]);
%% 
figure(2);
plot(t(1:n),z2s,'r',t(1:n),z222s,'b-',t(1:n),R11s,'m--',t(1:n),R12s,'m--','linewidth',2.5);xlabel('Time(sec)','FontName','Times New Roman');
legend('$e_2$ (The proposed nonlinear adaptive PI)','$e_2$ (Linear PI with time-varying gains by [15])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
xlabel('Time(sec)', 'FontName', 'Times New Roman');
ylabel(' $e_2$(rad)  ','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
  ylim([-10 10]);
  hold on
plot([T T], ylim, 'k--', 'LineWidth', 2);
legend('$e_2$ (The proposed nonlinear adaptive PI)','$e_2$ (Linear PI with time-varying gains by [15])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
axes('position',[0.55 0.15 0.35 0.25]);
plot(t(1:n),z2s,'r',t(1:n),z222s,'b-',t(1:n),R11s,'m--',t(1:n),R12s,'m--','linewidth',2.5);
xlim([2 20]);
ylim([-0.1 0.1]);
%% 
figure(3);
plot(t(1:n),u1s,'r',t(1:n),u111s,'b--','linewidth',2.5);xlabel('Time(sec)','FontName','Times New Roman');
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('$u_1$ (N$\cdot$m)','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
legend('$u_1$ (The proposed nonlinear adaptive PI)','$u_1$ (Linear PI with time-varying gains by [15])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
%% 
figure(4);
plot(t(1:n),u2s,'r',t(1:n),u222s,'b--','linewidth',2.5);xlabel('Time(sec)','FontName','Times New Roman');
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('$u_2$ ','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
legend('$u_2$ (The proposed nonlinear adaptive PI)','$u_2$ (Linear PI with time-varying gains by [15])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
figure(5);
plot(t(1:n),Q1s,'r',t(1:n),Y11s,'b--','linewidth',2.5);
hold on; % 保持当前图形
plot([t(1) t(n)], [-0.2 -0.2], 'k:', 'linewidth', 1.5); % 添加y=-0.2的水平线
plot([t(1) t(n)], [0.2 0.2], 'k:', 'linewidth', 1.5); % 添加y=0.2的水平线
hold off;
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('$u_I1$','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
legend('$\mathcal{N}_{I1}\big({\scriptstyle\int_0^tZ_1(\tau)d\tau},b_1\big)$ with constraint','$\int_0^tZ_1(\tau)d\tau$ without constraint','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
ylim([-0.5 0.5]);
%%
figure(6);
plot(t(1:n),Q2s,'r',t(1:n),Y21s,'b--','linewidth',2.5);
hold on; % 保持当前图形
plot([t(1) t(n)], [-0.2 -0.2], 'k:', 'linewidth', 1.5); % 添加y=-0.2的水平线
plot([t(1) t(n)], [0.2 0.2], 'k:', 'linewidth', 1.5); % 添加y=0.2的水平线
hold off;
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('$u_I2$','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
legend('$\mathcal{N}_{I2}\big({\scriptstyle\int_0^tZ_2(\tau)d\tau},b_2\big)$ with constraint','$\int_0^tZ_2(\tau)d\tau$ without constraint','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
ylim([-0.5 0.5]);
 % 
 figure(5);
plot(t(1:n),dz1s,'r',t(1:n),dz111s,'b-',t(1:n),R21s,'m--',t(1:n),R22s,'m--','linewidth',2.5);xlabel('Time(sec)','FontName','Times New Roman');
legend('$\dot e_1$ (The proposed nonlinear adaptive PI)','$\dot e_1$ (Linear PI with time-varying gains by [18])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
xlabel('Time(s)','FontName','Times New Roman','FontSize',12);
ylabel('$\dot e_1$(rad/s) ','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
 ylim([-5 5]);
axes('position',[0.55 0.15 0.35 0.25]);
plot(t(1:n),dz1s,'r',t(1:n),dz111s,'b-',t(1:n),R11s,'m--',t(1:n),R12s,'m--','linewidth',2.5);
xlim([4 20]);
ylim([-0.06 0.1]);
%% 
figure(6);
plot(t(1:n),dz2s,'r',t(1:n),dz222s,'b-',t(1:n),R21s,'m--',t(1:n),R22s,'m--','linewidth',2.5);xlabel('Time(sec)','FontName','Times New Roman');
legend('$\dot e_2$ (The proposed nonlinear adaptive PI)','$\dot e_2$ (Linear PI with time-varying gains by [18])','Interpreter','LaTex','FontName','Times New Roman','FontSize',10);
xlabel('Time(sec)', 'FontName', 'Times New Roman');
ylabel(' $\dot e_2$(rad/s)  ','Interpreter','LaTex','FontName','Times New Roman','FontSize',12);
 ylim([-5 5]);
axes('position',[0.55 0.15 0.35 0.25]);
plot(t(1:n),dz2s,'r',t(1:n),dz222s,'b-',t(1:n),R11s,'m--',t(1:n),R12s,'m--','linewidth',2.5);
xlim([4 20]);
ylim([-0.06 0.1]);