clear all
close all
warning off

T = 100;

n = 100;
% sigma = 1e-9;
sigma = 1e-4;
% sigma = 1e-3;

nTr1 = zeros(9,1);

nTr2 = zeros(9,1);

nTr3 = zeros(9,1);

nTr4 = zeros(9,1);

nTr5 = zeros(9,1);

for l=10:10:90

    mvecs(l/10) = l;
    g = linspace(1,n-l,n-l);


    %caso peggiore
    G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
    BestRkTr1(l/10) = sum(diag(G1(l+1:n,l+1:n)));


    %slow decay
    g2 = (-g+g(n-l)+1)./g(n-l);
    G2 = diag([ones(l,1); sigma*g2']);
    BestRkTr2(l/10) = sum(diag(G2(l+1:n,l+1:n)));


    %fast decay
    g3 = (1./(g).^5);
    G3 = diag([ones(l,1); sigma*g3']);
    BestRkTr3(l/10) = sum(diag(G3(l+1:n,l+1:n)));
    
    
    %slow decay
    g4 = (1./(g).^.5);
    G4 = diag([ones(l,1); sigma*g4']);
    BestRkTr4(l/10) = sum(diag(G4(l+1:n,l+1:n)));
   
    %caso migliore
    G5 = diag([ones(l,1); sigma; zeros(n-l-1,1)]);
    BestRkTr5(l/10) = sum(diag(G5(l+1:n,l+1:n)));

    for t=1:T
        
        [UBig,LhatBig] = PinvNystrom(G1,l);
        B1 = G1-UBig*LhatBig*UBig';
        nTr1(l/10) = nTr1(l/10) + trace(B1);
        
        [UBig,LhatBig] = PinvNystrom(G2,l);
        B2 = G2-UBig*LhatBig*UBig';
        nTr2(l/10) = nTr2(l/10) + trace(B2);
                
        [UBig,LhatBig] = PinvNystrom(G3,l);
        B3 = G3-UBig*LhatBig*UBig';
        nTr3(l/10) = nTr3(l/10) + trace(B3);

        [UBig,LhatBig] = PinvNystrom(G4,l);
        B4 = G4-UBig*LhatBig*UBig';
        nTr4(l/10) = nTr4(l/10) + trace(B4);

        [UBig,LhatBig] = PinvNystrom(G5,l);
        B5 = G5-UBig*LhatBig*UBig';
        nTr5(l/10) = nTr5(l/10) + trace(B5);

    end
    

end

nTr1 = 1/T * nTr1;

nTr2 = 1/T * nTr2;

nTr3 = 1/T * nTr3;

nTr4 = 1/T * nTr4;

nTr5 = 1/T * nTr5;

for tMV = 10:10:90

    g = linspace(1,n-tMV,n-tMV);

    %caso peggiore
    G1 = diag([ones(1,tMV), sigma*ones(1,n-tMV)]);
    BestRkTr1(tMV/10) = sum(diag(G1(tMV+1:n,tMV+1:n)));


    %slow decay
    g2 = (-g+g(n-tMV)+1)./g(n-tMV);
    G2 = diag([ones(tMV,1); sigma*g2']);
    BestRkTr2(tMV/10) = sum(diag(G2(tMV+1:n,tMV+1:n)));


    %fast decay
    g3 = (1./(g).^5);
    G3 = diag([ones(tMV,1); sigma*g3']);
    BestRkTr3(tMV/10) = sum(diag(G3(tMV+1:n,tMV+1:n)));
    
    
    %slow decay
    g4 = (1./(g).^.5);
    G4 = diag([ones(tMV,1); sigma*g4']);
    BestRkTr4(tMV/10) = sum(diag(G4(tMV+1:n,tMV+1:n)));
   
    %caso migliore
    G5 = diag([ones(tMV,1); sigma; zeros(n-tMV-1,1)]);
    BestRkTr5(tMV/10) = sum(diag(G5(tMV+1:n,tMV+1:n)));
end


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

l = 50;
g = linspace(1,n-l,n-l);


G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
    
g2 = (-g+g(n-l)+1)./g(n-l);
G2 = diag([ones(l,1); sigma*g2']);

g3 = (1./(g).^5);
G3 = diag([ones(l,1); sigma*g3']);

g4 = (1./(g).^.5);
G4 = diag([ones(l,1); sigma*g4']);
   
G5 = diag([ones(l,1); sigma; zeros(n-l-1,1)]);

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
