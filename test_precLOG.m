% 26 03 2025 Test Nystrom: come fare per evitare di fare Lanczos ma avere
% un risultato più accurato per shift piccoli.

clear all
close all

n = 100;
dl = 10;

Q = randn(n,n); [Q,~] = qr(Q);
g = linspace(1,n,n);
G = diag(1./(g).^2);
% G = diag(1./sqrt(g));
% G = diag(exp(-g));
A = Q*G*Q'; mi = 1e-3;
trA = sum(log(diag(G+mi*eye(n,n))),"all");

for l=10:10:10*dl
    [U,Lhat] = Nystrom(A,l); ll=Lhat(l,l);

    control1(l/10) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l)))) - trace(logm( (ll+mi)^(-1)*U*(Lhat+mi*eye(l))*U' + (eye(n) - U*U') ));
    P05 = ((ll+mi)*U*(Lhat+mi*eye(l))^-1*U' + (eye(n) - U*U') )^0.5;
    control2(l/10) = trace(logm(P05*(A+mi*eye(n))*P05));

    Est(l/10) = log((ll+mi)^(-l)*prod(diag(Lhat+mi*eye(l,l)))) + n*log(mi);
    
end

plot([10:10:10*dl],Est)
hold on
plot([10:10:10*dl],trA*ones(dl,1))