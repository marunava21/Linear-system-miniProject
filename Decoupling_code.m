clc
clear
A=[-1.7 -0.25 0; 23 -30 20; 0 -450 -740];
B=[5 0; -44 0; 0 -830];
C=[0 1 0;0 0 1];
D=[0 0; 0 0; 0 0];
x0=[1;100;200];

%intergrator - decoupled by state feedback
syms s;

sympref('FloatingPointOutput',true)
GS = C*inv(s*eye(3)-A)*B;

%Computing Bstar
J1 = C(1,:)*B;
J2 = C(2,:)*B ;
Bstar=[ J1(1,:) ; J2(1,:)];
CK1=det(Bstar); %check if Bstar is singular

%computing C star
I1= C(1,:)*A;
I2= C(2,:)*A;
Cstar = [I1(1,:) ; I2(1,:)];

%computing C double star with Pole 
I11= C(1,:)*A + 4*C(1,:)*eye(3);
I21= C(2,:)*A + 6*C(2,:)*eye(3);
Cdstar1 = ([I11(1,:) ; I21(1,:)]);

F = inv(Bstar);
K = F*Cstar;
K1 = (F*Cdstar1)

HS = C*inv(s*eye(3)-A + B*K)*B*F;

HS1 = C*inv(s*eye(3)-A + B*K1)*B*F;
simplify(HS1);
diag(HS1);
Acl = (A-B*K1);
E=eig(Acl)