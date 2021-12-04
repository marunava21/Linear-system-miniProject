clc
clear

""" statespace matrices  """;

A=[-1.7 -0.25 0; 23 -30 20; 0 -450 -740];
b=[5 0; -44 0; 0 -830];
c=[0 1 0;0 0 1];
D=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


"""   LQR   """ ;

Q=[0 0 0;0 0 0;0 0 10];
R=[1 0;0 1];
a1=A;
b1=-b*inv(R)*b';
c1=-Q;
d1=-A'; 
tau=[a1 b1;c1 d1];
[u,v]=eig(tau);
%%%%u has negative eigen values in 1,3 and 4th column, so we will divide 
%%the divide those eigen vectors into two parts which will be v and u of
%%the Pv=u.
v1=zeros(3,3);
u1=zeros(3,3);
count=1;
for i=1:6
    if(v(i,i)<0)
        v1(:,count)=real(u(1:3,i));
        u1(:,count)=real(u(4:6,i));
        count=count+1;
    end
end
P=u1*inv(v1);
K=inv(R)*b'*P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

""""""""""" observer design  """"""""""";

syms t1 t2 t3 g1 g2;
g1=10;
g2=10;
d=-20;
equ1=-1.7*t1+20*t2-d*t1;
equ2=-0.25*t1-30*t2-450*t3-d*t2-g1;
equ3=20*t2-740*t3-d*t3-g2;
[t1,t2,t3]=solve(equ1,equ2,equ3,t1,t2,t3);
t=[double(t1) double(t2) double(t3)];
e=t*b;
new_c1=[c;t];
obs_x=inv(new_c1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

""" Observer state and closed loop calculation """ ;

k11=double(K(1,1));k12=double(K(1,2));k13=double(K(1,3));
k21=double(K(2,1));k22=double(K(2,2));k23=double(K(2,3));
coef1=obs_x(1,1)*k11+k12;
coef2=obs_x(1,2)*k11+k13;
coef3=obs_x(1,3)*k11;
coef4=obs_x(1,1)*k21+k22;
coef5=obs_x(1,2)*k21+k23;
coef6=obs_x(1,3)*k21;
sympref('FloatingPointOutput',true);
syms y1 y2 sig;
u2=coef4*y1+coef5*y2+coef6*sig;
u1=coef1*y1+coef2*y2+coef3*sig;
sig_dot=d*sig+[1.3310 1.4110]*[u1;u2]+[1 1]*[y1;y2]
sig_dot_mat=coeffs(sig_dot)
u=[coef1 coef2 coef3;coef4 coef5 coef6];
u_new=b*u;
new_A=[-1.7 -0.25 0 0; 23 -30 20 0; 0 -450 -740 0;0 sig_dot_mat(1,2) sig_dot_mat(1,1) sig_dot_mat(1,3)];
Bu=[zeros(3,1) u_new;zeros(1,4)];
A_close=new_A+Bu;
new_c=[0 1 1 0];
eig(A_close)

