%%%TO BE FINISHED

function tr_est = TraceEst(A,l,t,m)

[n,~] = size(A);
B = (rand(n,t)<.5)*2 - 1; %Rademacher RM
[U,Lhat] = Nystrom(A,l);
x = zeros(n,t); tol=1e-8;
T = pcg_d(A,B,x,m,tol,U,Lhat);
[W,Ei]=eig(T);


