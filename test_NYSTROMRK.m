clear all

n = 100;
Q = randn(n,n); [Q,~] = qr(Q);
g=linspace(1,n,n);
G = diag(1./g.^0.5);
% G = diag(exp(-g));
A = Q*G*Q';
v=ones(n,1);
w=ones(n,1);
v=v/norm(v);
w=w/norm(w);
p.m=100;
p.tol=1e-10;
p.smin=G(end,end);
p.smax=G(1,1);
p.ch=0;
p.period=1;

% Bisogna capire quando si possono prendere poli non minuscoli...
% quando la matrice non ha decadimento degli autovalori...

X=lyap(-A,v*v');
[sz,Z]=Nystrom_RK_Lyap(-A,v,w,p,X);