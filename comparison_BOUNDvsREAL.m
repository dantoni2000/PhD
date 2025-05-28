% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off
T = 100;
n = 1000;
g = linspace(1,n,n)';
mi = 1;
ExperrtrBig=zeros(9,6);
Experrtr1=zeros(9,6);

for l=15:10:95

    for j = 1:6
        
        alpha = 0.64 +0.05*j;
        A = diag(alpha.^(g));
        trA(j) = sum(log(diag(A+mi*eye(n,n))),"all");

        mvecs((l-5)/10 ) = l;
        
        for t=1:T
        % Nystrom grande su A
        [UBig,LhatBig] = Nystrom(A,l);
        ExperrtrBig((l-5)/10,j) = ExperrtrBig((l-5)/10,j) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA(j));
        
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-3;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,3,1,1);
        Experrtr1((l-5)/10,j) = Experrtr1((l-5)/10,j) + abs(trr - trA(j));

        end
    
    end

end

etrBig = 1/T * ExperrtrBig;
etr1 = 1/T * Experrtr1;

for r=1:1:9
    figure(r)
    semilogy(0.69:0.05:0.94,etr1(r,:)','r')
    hold on
    semilogy(0.69:0.05:0.94,etrBig(r,:)','b')
    hold on
end

alpha = linspace(0.69,0.94,5000);

for tMV = 10:10:90

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;

        % boundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k+1) .* (sqrt(1-alpha.^(2*(n-k)))./sqrt(1-alpha.^2) + (k-1)/p .*sqrt(1-alpha.^(n-k))./ sqrt(1-alpha) + min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p), (k-1)/p) .*(1-alpha.^(n-k))./ (1-alpha));
        boundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k+3) .* (sqrt(1-alpha.^(2*(n-k-2)))./sqrt(1-alpha.^2) + (k+2-1)/p .*sqrt(1-alpha.^(n-k-2))./ sqrt(1-alpha) + min(exp(1) .* sqrt(k+2+p)/p .* sqrt((k+2-1)/p), (k+2-1)/p) .*(1-alpha.^(n-k-2))./ (1-alpha));
        boundBIGNys(:,s/2-1) = alpha.^(k+6) .* (1+(k+5-1)/(p)) .*(1-alpha.^(n-k-5))./(1-alpha); 
        badboundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k+3) .*(sqrt(1-alpha.^(2*(n-k-2)))./sqrt(1-alpha.^2) + (k+2-1)/p .*(1-alpha.^(n-k-2))./ (1-alpha));
    end
    
    BESTPrecSTE = min(boundPrecSTE');
    BESTBIGNys = min(boundBIGNys');
    BESTbadPrecSTE = min(badboundPrecSTE');
    
    tol = 1e-2;
    indexbound = abs(BESTPrecSTE - BESTBIGNys)./abs(BESTBIGNys) < tol;
    intersection = find(indexbound);
    
        if ~isempty(intersection) 
            coef = alpha(intersection(1));
            nu(tMV/10) = log(1/coef(1));
        else 
            coef = alpha(end);
            nu(tMV/10) = log(1/coef(1));
        end
        
        badindexbound = abs(BESTbadPrecSTE - BESTBIGNys)./abs(BESTBIGNys) < tol;
        badintersection = find(badindexbound);
    
        if ~isempty(badintersection) 
            badcoef = alpha(badintersection(1));
            badnu(tMV/10) = log(1/badcoef(1));
        else 
            badcoef = alpha(end);
            badnu(tMV/10) = log(1/badcoef(1));
        end
    
    figure(tMV/10)
    semilogy(alpha, BESTPrecSTE,'-or')
    hold on
    semilogy(alpha, BESTbadPrecSTE,'-dr')
    hold on
    semilogy(alpha, BESTBIGNys,'-ob')
    
end
figure(100)
    % plot(10:10:100,nu,'-*r')
    % hold on
    % plot(10:10:100,badnu,'-*b')
    % hold on
    plot(max(nu,badnu),'-*k')
    xlabel('$k$', 'Interpreter', 'latex')
    ylabel('$\nu$', 'Interpreter', 'latex')
    title('Coefficient $\nu$ for which $trace \log((I + x_{ii}^{-\nu}))$ can be better approximated with PrecSTE', 'Interpreter','latex')
