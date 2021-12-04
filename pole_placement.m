clc
clear
format short;
A=[-1.7 -0.25 0; 23 -30 20; 0 -450 -740];
b=[5 0; -44 0; 0 -830];
C=[0 1 0;0 0 1];
d=0;
x0=[1;100;200];
w=[b A*b A^2*b];
c=[w(:,1) w(:,3) w(:,2)];
c_inv=inv(c);
T=[c_inv(2,:);c_inv(2,:)*A;c_inv(3,:)];
T_inv=inv(T);
A_bar=T*A*T_inv;
B_bar=T*b;

%%%%%%% pole calculation %%%%%%%%%%%
sig=0.6;
w=0.25;
lambda=-(sig*w)+(w*(1-sig^2)^0.5)*i
syms s;
equ1=(s-(real(lambda)))^2+(imag(lambda))^2;
equ2=(s+5*real(-lambda));
pol_equ=expand(equ1*equ2);
pol_equ1=sym2poly(pol_equ);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Al=[0 1 0;0 0 1; -pol_equ1(1,4) -pol_equ1(1,3) -pol_equ1(1,2)];
kbar=mldivide(B_bar,(A_bar-Al));
% x0=[0;0;0];
k=kbar*T
sys_new=ss((A-b*k),b,C,d);
% initial(sys_new,x0);
[p,z]=pzmap(sys_new);
t=0:0.01:2;
u=[ones(size(t));zeros(size(t))];
[p,z]=pzmap(sys_new);
% step(sys_new);



