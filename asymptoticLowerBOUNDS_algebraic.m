% 08 04 2025 Test Hutchinson vs trace estimator: convergence of Lanczos and
% variance wrt samples

clear all
close all
warning off
T = 100;
n = 1000;
g = linspace(1,n,n)';
mi = 1;
ExperrtrBig=zeros(9,9);
Experrtr1=zeros(9,9);
ExperrtrNoPrec = zeros(9,9);

for l=15:10:95

    for j = 1:9
        
        alpha = .005 +j;
        A = diag(g.^-alpha);
        trA(j) = sum(log(diag(A+mi*eye(n,n))),"all");
        mvecs((l-5)/10 ) = l;
        
        for t=1:T
        % Nystrom grande su A
        [UBig,LhatBig] = Nystrom(A,l);
        ExperrtrBig((l-5)/10,j) = ExperrtrBig((l-5)/10,j) + abs(sum(log(diag(LhatBig+mi*eye(l,l))),"all") - trA(j))^2;
        
        %Nystrom con 1 Hutch 5 Lanczos
        l1 = l-5;
        [~,trr] = Nystrom_HUTCH(A,mi,l1,1,5,1,1);
        Experrtr1((l-5)/10,j) = Experrtr1((l-5)/10,j) + abs(trr - trA(j))^2;

        %Nystrom con n Hutch
        [~,trrNoPrec] = EST(A,mi,l/5,5,1);
        ExperrtrNoPrec((l-5)/10,j) = ExperrtrNoPrec((l-5)/10,j) + abs(trrNoPrec - trA(j))^2;

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

        % boundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k+1) .* (sqrt(1-alpha.^(2*(n-k)))./sqrt(1-alpha.^2) + (k-1)/p .*sqrt(1-alpha.^(n-k))./ sqrt(1-alpha) + min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p), (k-1)/p) .*(1-alpha.^(n-k))./ (1-alpha));        
        LOWERboundPrecSTE(:,s/2-1) = (1+ (1 + 1/4) .* min(1,((1+2*k/(p-1)).*(k+1).^(-alpha)+(2*exp(1)^2*(k+p)/(p^2-1).* 1./(alpha-1) .* ( (k.^(-alpha+1) - n.^(-alpha+1)) )))) ).*(1+1).^2;
        LOWERboundBIGNys(:,s/2-1) = ( (1+ 1/2 .* min(1,((1+2*(k+5)/(p-1)).*(k+6).^(-alpha)+(2*exp(1)^2*(k+5+p)/(p^2-1).* 1./(alpha-1) .* ( ((k+5).^(-alpha+1) - n.^(-alpha+1)) )))) ).*(1+1) ).^2; 

    end
    
    BESTLOWERboundPrecSTE = (min(LOWERboundPrecSTE')).^-1 .*2 .* 1./(2.*alpha-1) .* (tMV.^(-2*alpha+1) - n.^(-2*alpha+1));
    BESTLOWERBIGNys = (min(LOWERboundBIGNys')).^-1 .*(1./(alpha-1) .* ( (tMV+5).^(-alpha+1) - n.^(-alpha+1))).^2;
    BESTLOWERNONys = 10/(tMV+5) ./(1+1/2).^2 .* 1./(2.*alpha-1) .* (1.^(-2*alpha+1) - n.^(-2*alpha+1));

    tol = 1e-2;
    indexbound = abs(BESTLOWERboundPrecSTE - BESTLOWERBIGNys)./abs(BESTLOWERBIGNys) < tol;
    intersection = find(indexbound);
    
        if ~isempty(intersection) 
            coef = alpha(intersection(1));
            nu(tMV/10) = coef(1);
        else 
            coef = alpha(end);
            nu(tMV/10) = coef(1);
        end
        
    figure(tMV/10)
    semilogy(alpha, BESTLOWERboundPrecSTE,'-or')
    hold on
    semilogy(alpha, BESTLOWERBIGNys,'-*c')
    hold on
    semilogy(alpha, BESTLOWERNONys,'-dg')

    xlabel('$\alpha$', 'Interpreter','latex')
    ylabel('error')
    legend('error Nystrom+Lanczos', 'error Big Nystrom', 'error No Nystrom', 'Lower Bound Nystrom + Lanczos', 'Lower Bound Big Nystrom', 'Lower Bound NO Nystrom')
    
end
figure(100)
    % plot(10:10:100,nu,'-*r')
    % hold on
    % plot(10:10:100,badnu,'-*b')
    % hold on
    plot(nu,'-*k')
    xlabel('$k$', 'Interpreter', 'latex')
    ylabel('$\alpha$', 'Interpreter', 'latex')
    title('Coefficient $\log(\alpha)$ for which $trace \log((I + \alpha^{-x_{ii}}))$ can be better approximated with PrecSTE', 'Interpreter','latex')
