clc
clear
%""" The augmented system """%
A=[-1.7 -0.25 0 0 0;23 -30 20 0 0;0 -450 -740 0 0;0 -1 0 0 0;0 0 -1 0 0];
b=[5 0;-44 0;0 -830;0 0;0 0];
c=[0 1 0 0 0;0 0 1 0 0];
d=0;
bw=[-2;5;0;0;0];
br=[0;0;0;0;1];


%""" LQR """%
Q=[1 0 0 0 0;0 5 0 0 0;0 0 5 0 0;0 0 0 5 0;0 0 0 0 5];
R=[1 0;0 1];
a1=A;
b1=-b*inv(R)*b';
c1=-Q;
d1=-A'; 
tau=[a1 b1;c1 d1];
[u,v]=eig(tau);
v1=zeros(5,5);
u1=zeros(5,5);
count=1;
for i=1:10
    if(v(i,i)<0)
        v1(:,count)=real(u(1:5,i));
        u1(:,count)=real(u(6:10,i));
        count=count+1;
    end
end
P=u1*inv(v1);
K=inv(R)*b'*P