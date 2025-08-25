clear all
close all
warning off

T = 50;

n = 100;
% sigma = 1e-9;
sigma = 1e-6;
% sigma = 1e-3;

nTr1 = zeros(9,1);
nalmostTr1 = zeros(9,1);

nTr2 = zeros(9,1);
nalmostTr2 = zeros(9,1);

nTr3 = zeros(9,1);
nalmostTr3 = zeros(9,1);

nTr4 = zeros(9,1);
nalmostTr4 = zeros(9,1);

nTr5 = zeros(9,1);
nalmostTr5 = zeros(9,1);

for l=10:10:90

    mvecs(l/10) = l;
    g = linspace(1,l,l)';
    almostg = linspace(1,l-1,l-1)';


    %caso peggiore
    G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
    almostG1 = diag([ones(1,l-1), sigma*ones(1,n-l+1)]);
    
    BestRkTr1(l/10) = sum(diag(G1(l+1:n,l+1:n)));
    BestRkalmostTr1(l/10) = sum(diag(almostG1(l+1:n,l+1:n)));

 
    %slow decay
    g2 = (-g'+g(l,1)+1)./g(l,1);
    G2 = diag([max(g2,sigma), sigma*ones(1,n-l)]);
    almostg2 = (-almostg'+almostg(l-1,1)+1)./almostg(l-1,1);
    almostG2 = diag([max(almostg2,sigma), sigma*ones(1,n-l+1)]);

    BestRkTr2(l/10) = sum(diag(G2(l+1:n,l+1:n)));
    BestRkalmostTr2(l/10) = sum(diag(almostG2(l+1:n,l+1:n)));


    %fast decay
    g3 = (1./(g').^5);
    G3 = diag([max(g3,sigma), sigma*ones(1,n-l)]);
    almostg3 = (1./(almostg').^5);
    almostG3 = diag([max(almostg3,sigma), sigma*ones(1,n-l+1)]);

    BestRkTr3(l/10) = sum(diag(G3(l+1:n,l+1:n)));
    BestRkalmostTr3(l/10) = sum(diag(almostG3(l+1:n,l+1:n)));


    %slow decay
    g4 = (1./(g').^.5);
    G4 = diag([max(g4,sigma), sigma*ones(1,n-l)]);
    almostg4 = (1./(almostg').^.5);
    almostG4 = diag([max(almostg4,sigma), sigma*ones(1,n-l+1)]);

    BestRkTr4(l/10) = sum(diag(G4(l+1:n,l+1:n)));
    BestRkalmostTr4(l/10) = sum(diag(almostG4(l+1:n,l+1:n)));


    %caso migliore
    G5 = diag([1, sigma*ones(1,n-1)]);
    
    BestRkTr5(l/10) = sum(diag(G5(l+1:n,l+1:n)));

    for t=1:T
        
        [UBig,LhatBig] = PinvNystrom(G1,l);
        B1 = G1-UBig*LhatBig*UBig';
        nTr1(l/10) = nTr1(l/10) + trace(B1);
        
        [UBig,LhatBig] = PinvNystrom(almostG1,l);
        almostB1 = almostG1-UBig*LhatBig*UBig';
        nalmostTr1(l/10) = nalmostTr1(l/10) + trace(almostB1);

        [UBig,LhatBig] = PinvNystrom(G2,l);
        B2 = G2-UBig*LhatBig*UBig';
        nTr2(l/10) = nTr2(l/10) + trace(B2);
                
        [UBig,LhatBig] = PinvNystrom(almostG2,l);
        almostB2 = almostG2-UBig*LhatBig*UBig';
        nalmostTr2(l/10) = nalmostTr2(l/10) + trace(almostB2);
    
        [UBig,LhatBig] = PinvNystrom(G3,l);
        B3 = G3-UBig*LhatBig*UBig';
        nTr3(l/10) = nTr3(l/10) + trace(B3);

        [UBig,LhatBig] = PinvNystrom(almostG3,l);
        almostB3 = almostG3-UBig*LhatBig*UBig';
        nalmostTr3(l/10) = nalmostTr3(l/10) + trace(almostB3);
    
        [UBig,LhatBig] = PinvNystrom(G4,l);
        B4 = G4-UBig*LhatBig*UBig';
        nTr4(l/10) = nTr4(l/10) + trace(B4);

        [UBig,LhatBig] = PinvNystrom(almostG4,l);
        almostB4 = almostG4-UBig*LhatBig*UBig';
        nalmostTr4(l/10) = nalmostTr4(l/10) + trace(almostB4);

        [UBig,LhatBig] = PinvNystrom(G5,l);
        B5 = G5-UBig*LhatBig*UBig';
        nTr5(l/10) = nTr5(l/10) + trace(B5);

    end
    

end

nTr1 = 1/T * nTr1;
nalmostTr1 = 1/T * nalmostTr1;

nTr2 = 1/T * nTr2;
nalmostTr2 = 1/T * nalmostTr2;

nTr3 = 1/T * nTr3;
nalmostTr3 = 1/T * nalmostTr3;

nTr4 = 1/T * nTr4;
nalmostTr4 = 1/T * nalmostTr4;

nTr5 = 1/T * nTr5;

for tMV = 10:10:90

    g = linspace(1,tMV,tMV)';
    almostg = linspace(1,tMV-1,tMV-1)';
    G1 = diag([ones(1,tMV), sigma*ones(1,n-tMV)]);
    almostG1 = diag([ones(1,tMV-1), sigma*ones(1,n-tMV+1)]);
    
    BestRkTr1(tMV/10) = sum(diag(G1(tMV+1:n,tMV+1:n)));
    BestRkalmostTr1(tMV/10) = sum(diag(almostG1(tMV+1:n,tMV+1:n)));


    g2 = (-g'+g(tMV,1)+1)./g(tMV,1);
    G2 = diag([max(g2,sigma), sigma*ones(1,n-tMV)]);
    almostg2 = (-almostg'+almostg(tMV-1,1)+1)./almostg(tMV-1,1);
    almostG2 = diag([max(almostg2,sigma), sigma*ones(1,n-tMV+1)]);

    BestRkTr2(tMV/10) = sum(diag(G2(tMV+1:n,tMV+1:n)));
    BestRkalmostTr2(tMV/10) = sum(diag(almostG2(tMV+1:n,tMV+1:n)));


    g3 = (1./(g').^5);
    G3 = diag([max(g3,sigma), sigma*ones(1,n-tMV)]);
    almostg3 = (1./(almostg').^.5);
    almostG3 = diag([max(almostg3,sigma), sigma*ones(1,n-tMV+1)]);

    BestRkTr3(tMV/10) = sum(diag(G3(tMV+1:n,tMV+1:n)));
    BestRkalmostTr3(tMV/10) = sum(diag(almostG3(tMV+1:n,tMV+1:n)));


    g4 = (1./(g').^.5);
    G4 = diag([max(g4,sigma), sigma*ones(1,n-tMV)]);
    almostg4 = (1./(almostg').^.5);
    almostG4 = diag([max(almostg4,sigma), sigma*ones(1,n-tMV+1)]);

    BestRkTr4(tMV/10) = sum(diag(G4(tMV+1:n,tMV+1:n)));
    BestRkalmostTr4(tMV/10) = sum(diag(almostG4(tMV+1:n,tMV+1:n)));

    G5 = diag([1, sigma*ones(1,n-1)]);
    BestRkTr5(tMV/10) = sum(diag(G5(tMV+1:n,tMV+1:n)));
end


figure(2)
semilogy(mvecs,nTr1'./BestRkTr1,'-db')
hold on
semilogy(mvecs,nalmostTr1'./BestRkalmostTr1,'-ob')
hold on
semilogy(mvecs,nTr2'./BestRkTr2,'-dc')
hold on
semilogy(mvecs,nalmostTr2'./BestRkalmostTr2,'-oc')
hold on
semilogy(mvecs,nTr3'./BestRkTr3,'-dm')
hold on
semilogy(mvecs,nalmostTr3'./BestRkalmostTr3,'-om')
hold on
semilogy(mvecs,nTr4'./BestRkTr4,'-dr')
hold on
semilogy(mvecs,nalmostTr4'./BestRkalmostTr4,'-or')
hold on
semilogy(mvecs,nTr5'./BestRkTr5,'-dk')

xlabel('MatVecs')
ylabel('error')
title('Comparison of the bounds for Nystrom in Trace Norm')
legend('case 1', 'case almost1', 'case 2', 'case almost 2', 'case 3', 'case almost 3', 'case 4', 'case almost 4', 'case 5', 'interpreter', 'Latex')

l = 50;
g = linspace(1,l,l)';
almostg = linspace(1,l-1,l-1)';
G1 = diag([ones(1,l), sigma*ones(1,n-l)]);
almostG1 = diag([ones(1,l-1), sigma*ones(1,n-l+1)]);


g2 = (-g'+g(l,1)+1)./g(l,1);
G2 = diag([max(g2,sigma), sigma*ones(1,n-l)]);
almostg2 = (-almostg'+almostg(l-1,1)+1)./almostg(l-1,1);
almostG2 = diag([max(almostg2,sigma), sigma*ones(1,n-l+1)]);


g3 = (1./(g').^5);
G3 = diag([max(g3,sigma), sigma*ones(1,n-l)]);
almostg3 = (1./(almostg').^5);
almostG3 = diag([max(almostg3,sigma), sigma*ones(1,n-l+1)]);


g4 = (1./(g').^.5);
G4 = diag([max(g4,sigma), sigma*ones(1,n-l)]);
almostg4 = (1./(almostg').^.5);
almostG4 = diag([max(almostg4,sigma), sigma*ones(1,n-l+1)]);


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
