% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off
T = 100;
n = 2000;
g = linspace(1,n,n)';
mi = 1;
ExperrtrBig=zeros(19,1);
Experrtr1=zeros(19,1);

A = diag(sqrt(rand(n,1).^.2)); A = diag(sort(diag(A),'descend')); A = A./norm(A);
trA = sum(log(diag(A+mi*eye(n,n))),"all");

for l=15:10:195

    mvecs((l-5)/10 ) = l;
            
    for t=1:T
    % Nystrom grande su A
    [UBig,LhatBig] = Nystrom(A,l);
    ExperrtrBig((l-5)/10) = ExperrtrBig((l-5)/10) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA);
        
    %Nystrom con 1 Hutch 5 Lanczos
    l1 = l-5;
    [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
    Experrtr1((l-5)/10) = Experrtr1((l-5)/10) + abs(trr - trA);
            
    end
    

end

etrBig = 1/T * ExperrtrBig;
etr1 = 1/T * Experrtr1;

for tMV = 10:10:190

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;
        
        boundBIGNys(:,s/2-1) = ((1 + (k+5)/(p-1)).* norm(A(k+1:n,k+1:n).^0.5 ,'fro')^2);
        badboundPrecSTE(:,s/2-1) = sqrt(2) .* ((1 + ((k)/(p-1))).* norm(A(k+1:n,k+1:n),'fro') );
    end
    
    BESTBIGNys(tMV/10) = min(boundBIGNys');
    BESTbadPrecSTE(tMV/10) = min(badboundPrecSTE');
    
    
    
end

figure(1)
semilogy(mvecs,etr1,'r')
hold on
semilogy(mvecs,etrBig,'b')
hold on
semilogy(mvecs,BESTbadPrecSTE','-dr')
hold on
semilogy(mvecs,BESTBIGNys','-ob')

xlabel('MatVecs')
ylabel('error')
legend('error Nystrom+Lanczos', 'error Big Nystrom', 'conjecture Nystrom+Lanczos', 'bound Big Nystrom')

figure(100)
plot(diag(A))