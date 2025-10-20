clear all
close all
warning off

T = 1000;

n = 2000;
% sigma = 1e-9;
sigma = 1e-7;
% sigma = 1e-3;

nTr1 = zeros(1,T);

l = 1000;

G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
BestRkTr1 = sum(diag(G1(l+1:n,l+1:n)));
denom = 1/l + 3 * sigma *(n-l);


for t=1:T
    [UBig,LhatBig] = PinvNystrom(G1,l);
    B1 = G1-UBig*LhatBig*UBig';
    nTr1(1,t) = trace(B1);
end

figure(1)
subplot('Position', [0.05 0.3 0.4 0.5])
semilogy(1:T,nTr1'./BestRkTr1,'-db')
hold on
semilogy(1:T,ones(T,1)./denom,'-dc')
xlabel('$t$','interpreter','Latex')
ylabel('ratio')
title('Comparison of the bounds for Nystrom in Trace Norm')
legend('computed ratio', 'lower bound','interpreter','Latex')

figure(1)
subplot('Position', [0.55 0.3 0.4 0.5])
semilogy(1:T,(nTr1'./BestRkTr1).*denom,'-dk')
hold on
semilogy(1:T, ones(T,1).*sum((nTr1'./BestRkTr1).*denom)/T, 'r')
xlabel('try')
ylabel('error')
title('ratio between the computed quantity and the lower bound')
legend('ratio', 'interpreter', 'Latex')

sum(((nTr1'./BestRkTr1).*denom<1))