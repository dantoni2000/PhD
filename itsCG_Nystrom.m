clear all, close all

c = 5;

ee = 1e-8;
n = 1000;

P = 1:50;

ord = 0.5;
Q0 = randn(n,n); [Q0,~] = qr(Q0);
g = linspace(1,n,n)';
mi = 1e-4;
G0 = diag(1./(g).^ord);
A0 = Q0*G0*Q0';
b0 = A0*ones(n,1);
kA0 = (G0(1,1)+mi)/(G0(end,end)+mi);
K0 = log(2*kA0/ee);

for s=1:50
    [U0,Lhat0] = Nystrom(A0,2*s-1);
    P0s = U0*(Lhat0+mi*eye(2*s-1))^0.5*U0' + (mi)^0.5*(eye(n) - U0*U0');
    [~,~,~,its0(s)] = pcg(A0 + mi*eye(n),b0,1e-10,n,P0s,P0s);
end

CG_ord = K0*sqrt(0.5/mi)*sqrt(mi + 4 / (1-ord) * (n^(1-ord) - P.^(1-ord)));
figure(1)
plot (2*P-1, 2*P-1 + CG_ord,'g')
hold on
plot (2*P-1, 2*P-1 + its0,'o-g')


Q = randn(n,n); [Q,~] = qr(Q);
mi = 1e-4;
G = diag(1./(g));
A = Q*G*Q';
b = A*ones(n,1);
kA = (G(1,1)+mi)/(G(end,end)+mi);
K1 = log(2*kA/ee);

for s=1:50
    [U,Lhat] = Nystrom(A,2*s-1);
    Ps = U*(Lhat+mi*eye(2*s-1))^0.5*U' + (mi)^0.5*(eye(n) - U*U');
    [~,~,~,its(s)] = pcg(A + mi*eye(n),b,1e-10,n,Ps,Ps);
end

CG_l = K1*sqrt(0.5/mi)*sqrt(mi + 4*log(n./P));
figure(2)
plot (2*P-1, 2*P-1 + CG_l,'r')
hold on
plot (2*P-1, 2*P-1 + its,'o-r')

Q2 = randn(n,n); [Q2,~] = qr(Q2);
G2 = diag(1./(g.^2));
mi = 1;
A2 = Q2*G2*Q2';
b2 = A2*ones(n,1);
kA2 = (G2(1,1)+mi)/(G2(end,end)+mi);
K2 = log(2*kA2/ee);

for s=1:50
    [U2,Lhat2] = Nystrom(A2,2*s-1);
    P2s = U2*(Lhat2+mi*eye(2*s-1))^0.5*U2' + (mi)^0.5*(eye(n) - U2*U2');
    [~,~,~,its2(s)] = pcg(A2 + mi*eye(n),b2,1e-10,n,P2s,P2s);
end

CG_q = K2*sqrt(0.5/mi)*sqrt(mi + 4./P - 4/n);
figure(3)
plot (2*P-1, 2*P-1 + CG_q,'b')
hold on
plot (2*P-1, 2*P-1 + its2,'o-b')

Q3 = randn(n,n); [Q3,~] = qr(Q3);
G3 = diag(1./(g.^3));
mi = 1e-4;
A3 = Q3*G3*Q3';
b3 = A3*ones(n,1);
kA3 = (G3(1,1)+mi)/(G3(end,end)+mi);
K3 = log(2*kA3/ee);

for s=1:50
    [U3,Lhat3] = Nystrom(A3,2*s-1);
    P3s = U3*(Lhat3+mi*eye(2*s-1))^0.5*U3' + (mi)^0.5*(eye(n) - U3*U3');
    [~,~,~,its3(s)] = pcg(A3 + mi*eye(n),b3,1e-10,n,P3s,P3s);
end

CG_c = K3*sqrt(0.5/mi)*sqrt(mi + 4/2 * (1./P.^2 - 1/n^2) + 4/6* (1./(P.^3)-1./(n^3)));
figure(4)
plot (2*P-1, 2*P-1 + CG_c,'k')
hold on
plot (2*P-1, 2*P-1 + its3,'o-k')

Q4 = randn(n,n); [Q4,~] = qr(Q4);
G4 = diag(c.^(-g));
mi = 1e-4;
A4 = Q4*G4*Q4';
b4 = A4*ones(n,1);
kA4 = (G4(1,1)+mi)/(G4(end,end)+mi);
K4 = log(2*kA4/ee);

for s=1:50
    [U4,Lhat4] = Nystrom(A4,2*s-1);
    P4s = U4*(Lhat4+mi*eye(2*s-1))^0.5*U4' + (mi)^0.5*(eye(n) - U4*U4');
    [~,~,~,its4(s)] = pcg(A4 + mi*eye(n),b4,1e-10,n,P4s,P4s);
end

CG_exp = K4*sqrt(0.5/mi)*(mi+4*(c.^-(P)).^0.5/sqrt(1-1/c));
figure(5)
plot(2*P-1, 2*P-1 + CG_exp,'m')
hold on
plot(2*P-1, 2*P-1 + its4,'o-m')

b5 = randn(n,1);
kA4 = (G4(1,1)+mi)/(G4(end,end)+mi);
K4 = log(2*kA4/ee);

for s=1:50
    [U4,Lhat4] = Nystrom(A4,2*s-1);
    P4s = U4*(Lhat4+mi*eye(2*s-1))^0.5*U4' + (mi)^0.5*(eye(n) - U4*U4');
    [~,~,~,its5(s)] = pcg(A4 + mi*eye(n),b5,1e-10,n,P4s,P4s);
end

CG_exp = K4*sqrt(0.5/mi)*(mi+4*(c.^-(P)).^0.5/sqrt(1-1/c));
figure(6)
plot(2*P-1, 2*P-1 + CG_exp,'m')
hold on
plot(2*P-1, 2*P-1 + its5,'o-m')

% P = 1:350;
% 
% n = 700;
% Q = randn(n,n); [Q,~] = qr(Q);
% g = linspace(1,n,n)';
% mi = 1;
% G = 10^5*diag((-g+g(n,1)+1)./g(n,1)).^10;
% A = Q*G*Q';