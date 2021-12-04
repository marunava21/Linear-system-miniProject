clc
clear
A=[-1.7 -0.25 0; 23 -30 20; 0 -450 -740];
b=[5 0; -44 0; 0 -830];
C=[0 1 0;0 0 1];
d=0;
Q=[0 0 0;0 60 0;0 0 900];
R=[10 0;0 10];
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
sys=ss((A-b*K),b,C,d);
step(sys);
l=stepinfo(sys);
l(1,1)
l(2,2)
