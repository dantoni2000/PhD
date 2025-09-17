clear all
close all
warning off

T = 10;

n = 1000;
% sigma = 1e-9;
sigma = 1e-6;
% sigma = 1e-3;

nTr1 = zeros(49,T);

nTr2 = zeros(49,T);

nTr3 = zeros(49,T);

nTr4 = zeros(49,T);

nTr5 = zeros(49,T);

for l=10:10:490

    mvecs(l/10) = l;
    g = linspace(1,l,l)';

    %caso peggiore
    G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
    BestRkTr1(l/10) = sum(diag(G1(l+1:n,l+1:n)));


    %slow decay
    g2 = (-g'+g(l,1)+1)./g(l,1);
    G2 = diag([max(g2,sigma), sigma*ones(1,n-l)]);
    % G2 = diag([max(g2,sigma)/(sum(max(g2,sigma)))*l, sigma*ones(1,n-l)]);
    BestRkTr2(l/10) = sum(diag(G2(l+1:n,l+1:n)));


    %medium decay
    % g3 = (1./(g').^5);
    g3 = (1./(g').^.5);
    G3 = diag([ones(1,l/2), ones(1,l/2)*g3(1,l), sigma*ones(1,n-l)]);
    % G3 = diag([max(g3,sigma)/(sum(max(g3,sigma)))*l, sigma*ones(1,n-l)]);
    BestRkTr3(l/10) = sum(diag(G3(l+1:n,l+1:n)));


    %slow decay
    g4 = (1./(g').^.5);
    G4 = diag([max(g4,sigma), sigma*ones(1,n-l)]);
    % G4 = diag([max(g4,sigma)/(sum(max(g4,sigma)))*l, sigma*ones(1,n-l)]);
    BestRkTr4(l/10) = sum(diag(G4(l+1:n,l+1:n)));


    %caso Wang Zhang
    sigma2 = 1 - sigma;
    % sigma2 = sigma;
    G5 = sigma2 * ones(n,n); G5 = G5 + eye(n,n)*(1-sigma2);
    % G5 = diag([l - (l-1)*sigma, sigma*ones(1,n-1)]);    
    BestRkTr5(l/10) = (n-l)*(1-sigma2);

    %caso Derezinski
    G5 = construct_simplex_matrix(l, 1e-4);
    L5 = svd(G5);
    BestRkTr5(l/10) = L5(end);
    
    % %caso migliore
    % G5 = diag([1, sigma*ones(1,n-1)]);
    % % G5 = diag([l - (l-1)*sigma, sigma*ones(1,n-1)]);    
    % BestRkTr5(l/10) = sum(diag(G5(l+1:n,l+1:n)));

    for t=1:T
        
        [UBig,LhatBig] = PinvNystrom(G1,l);
        B1 = G1-UBig*LhatBig*UBig';
        nTr1(l/10,t) = trace(B1);
        
        [UBig,LhatBig] = PinvNystrom(G2,l);
        B2 = G2-UBig*LhatBig*UBig';
        nTr2(l/10,t) = trace(B2);         
    
        [UBig,LhatBig] = PinvNystrom(G3,l);
        B3 = G3-UBig*LhatBig*UBig';
        nTr3(l/10,t) = trace(B3);

        [UBig,LhatBig] = PinvNystrom(G4,l);
        B4 = G4-UBig*LhatBig*UBig';
        nTr4(l/10,t) = trace(B4);

        [UBig,LhatBig] = PinvNystrom(G5,l);
        B5 = G5-UBig*LhatBig*UBig';
        nTr5(l/10,t) = trace(B5);

    end
    

end

nTr1 = median(nTr1,2);
% nTr1 = trimmean(nTr1,10,2);

nTr2 = median(nTr2,2);
% nTr2 = trimmean(nTr2,10,2);

nTr3 = median(nTr3,2);
% nTr3 = trimmean(nTr3,10,2);

nTr4 = median(nTr4,2);
% nTr4 = trimmean(nTr4,10,2);

nTr5 = median(nTr5,2);
% nTr5 = trimmean(nTr5,10,2);

% for tMV = 10:10:490
% 
%     g = linspace(1,tMV,tMV)';
%     G1 = diag([ones(1,tMV-5), sigma*ones(1,n-tMV+5)]);
%     BestRkTr1(tMV/10) = sum(diag(G1(tMV+1:n,tMV+1:n)));
% 
%     g2 = (-g'+g(tMV,1)+1)./g(tMV,1);
%     G2 = diag([max(g2,sigma), sigma*ones(1,n-tMV)]);
%     % G2 = diag([max(g2,sigma)/(sum(max(g2,sigma)))*tMV, sigma*ones(1,n-tMV)]);
%     BestRkTr2(tMV/10) = sum(diag(G2(tMV+1:n,tMV+1:n)));
% 
%     % g3 = (1./(g').^5);
%     g3 = (1./(g').^.5);
%     G3 = diag([ones(1,tMV/2), ones(1,tMV/2)*g3(1,tMV), sigma*ones(1,n-tMV)]);
%     % G3 = diag([max(g3,sigma)/(sum(max(g3,sigma)))*tMV, sigma*ones(1,n-tMV)]);
%     BestRkTr3(tMV/10) = sum(diag(G3(tMV+1:n,tMV+1:n)));
% 
%     g4 = (1./(g').^.5);
%     G4 = diag([max(g4,sigma), sigma*ones(1,n-tMV)]); 
%     % G4 = diag([max(g4,sigma)/(sum(max(g4,sigma)))*tMV, sigma*ones(1,n-tMV)]); 
%     BestRkTr4(tMV/10) = sum(diag(G4(tMV+1:n,tMV+1:n)));
% 
%     G5 = diag([1, sigma*ones(1,n-1)]);
%     % G5 = diag([tMV - (tMV-1)*sigma, sigma*ones(1,n-1)]);
%     BestRkTr5(tMV/10) = sum(diag(G5(tMV+1:n,tMV+1:n)));
% end


figure(2)
semilogy(mvecs,nTr1'./BestRkTr1,'-db')
hold on
semilogy(mvecs,nTr2'./BestRkTr2,'-dc')
hold on
semilogy(mvecs,nTr3'./BestRkTr3,'-dm')
hold on
semilogy(mvecs,nTr4'./BestRkTr4,'-dr')
hold on
semilogy(mvecs,nTr5'./BestRkTr5,'-dk')

xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Trace Norm')
legend('case 1', 'case 2', 'case 3', 'case 4', 'case 5', 'interpreter', 'Latex')

l = 490;
g = linspace(1,l,l)';
G1 = diag([ones(1,l), sigma*ones(1,n-l)]);

g2 = (-g'+g(l,1)+1)./g(l,1);
G2 = diag([max(g2,sigma), sigma*ones(1,n-l)]);

g4 = (1./(g').^.5);
G3 = diag([ones(1,l/2),ones(1,l/2)*g4(1,l), sigma*ones(1,n-l)]);

g4 = (1./(g').^.5);
G4 = diag([max(g4,sigma), sigma*ones(1,n-l)]);

G5 = diag([1, sigma*ones(1,n-1)]);

figure
semilogy(diag(G1),'b')
hold on
semilogy(diag(G2),'c')
hold on
semilogy(diag(G3),'m')
hold on
semilogy(diag(G4),'r')
hold on
semilogy(diag(G5),'k')
legend('case 1', 'case 2', 'case 3', 'case 4', 'case 5', 'interpreter', 'Latex')
