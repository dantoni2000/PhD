% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off
T = 100;
n = 1000;
g = linspace(1,n,n)';
mi = 1;
ExperrtrBig = zeros(9,9);
Experrtr1 = zeros(9,9);
ExperrtrNoPrec = zeros(9,9);

for l=15:10:95

    for j = 1:9
        alpha = .5 +j;
        A = diag(g.^-alpha);
        trA(j) = sum(log(diag(A+mi*eye(n,n))),"all");

        % norm(A(l+1:n,l+1:n),'fro') - 1./sqrt(2.*alpha-1) .* sqrt( (l.^(-2*alpha+1) - n.^(-2*alpha+1)) )
        % norm(A(l+1:n,l+1:n).^0.5 ,'fro')^2 - 1./(alpha-1) .* ( (l.^(-alpha+1) - n.^(-alpha+1)) )
        % norm(A(l+1:n,l+1:n)) - (l+1).^ (-alpha)

        mvecs((l-5)/10 ) = l;
        
        for t=1:T
        % Nystrom grande su A
        [UBig,LhatBig] = Nystrom(A,l);
        ExperrtrBig((l-5)/10,j) = ExperrtrBig((l-5)/10,j) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA(j));
    
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
        Experrtr1((l-5)/10,j) = Experrtr1((l-5)/10,j) + abs(trr - trA(j));
        
        %Nystrom con n Hutch
        [~,trrNoPrec] = EST(A,mi,l/5,5,1);
        ExperrtrNoPrec((l-5)/10,j) = ExperrtrNoPrec((l-5)/10,j) + abs(trrNoPrec - trA(j));

        end
    
    end

end

etrBig = 1/T * ExperrtrBig;
etr1 = 1/T * Experrtr1;
etrNoPrec = 1/T * ExperrtrNoPrec;

for r=1:1:9
    figure(r)
    semilogy(1.5:1:9.5,etr1(r,:)','r')
    hold on
    semilogy(1.5:1:9.5,etrBig(r,:)','b')
    hold on
    semilogy(1.5:1:9.5,etrNoPrec(r,:)','g')
    hold on
end

alpha = linspace(1.5,9.5,5000);
for tMV = 10:10:90

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;
        
        boundPrecSTE(:,s/2-1) = sqrt(2) .* (1./sqrt(2.*alpha-1) .* sqrt( (k.^(-2*alpha+1) - n.^(-2*alpha+1)) ) + (k)/(p-1) ./sqrt(alpha-1).*(k+1).^(-alpha/2).*sqrt( (k.^(-alpha+1) - n.^(-alpha+1)) ) + (k.^(-alpha+1) - n.^(-alpha+1)) .* min(exp(1) .* sqrt(k+p)/p .* sqrt((k)/(p-1)),(k)/(p-1)) .* 1./(alpha-1));
        boundBIGNys(:,s/2-1) = ((1 + (k+5)/(p-1)).*1./(alpha-1).*((k+5).^(-alpha+1) - n.^(-alpha+1)));
        badboundPrecSTE(:,s/2-1) = sqrt(2) .* ((1 + ((k)/(p-1))).*1./sqrt(2.*alpha-1).* sqrt( (k.^(-2*alpha+1) - n.^(-2*alpha+1)) ));
    end
    
    BESTPrecSTE = min(boundPrecSTE');
    BESTBIGNys = min(boundBIGNys');
    BESTbadPrecSTE = min(badboundPrecSTE');
    BESTNONys = sqrt( 10/(tMV+5) .* ( zeta(2.*alpha) - 1./(2.*alpha-1) .* n.^(-2*alpha+1)) );
    
    tol = 1e-3;
    indexbound = abs(BESTPrecSTE - BESTBIGNys)./abs(BESTBIGNys) < tol;
    intersection = find(indexbound);
    
        if ~isempty(intersection) 
            coef = alpha(intersection(1));
            nu(tMV/10) = coef(1);
        else 
            coef = alpha(end);
            nu(tMV/10) = coef(1);
        end
        
        badindexbound = abs(BESTbadPrecSTE - BESTBIGNys)./abs(BESTBIGNys) < tol;
        badintersection = find(badindexbound);
    
        if ~isempty(badintersection) 
            badcoef = alpha(badintersection(1));
            badnu(tMV/10) = badcoef(1);
        else 
            badcoef = alpha(end);
            badnu(tMV/10) = badcoef(1);
        end
    
    figure(tMV/10)
    semilogy(alpha, BESTNONys,'-og')
    hold on
    semilogy(alpha, BESTbadPrecSTE,'-dr')
    hold on
    semilogy(alpha, BESTBIGNys,'-ob')

    xlabel('$\nu$', 'Interpreter','latex')
    ylabel('error')
    legend('error Nystrom+Lanczos', 'error Big Nystrom','bound Nystrom+Lanczos', 'conjecture Nystrom+Lanczos', 'bound Big Nystrom')
    
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
