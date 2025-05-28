clear, close all

n = 1000;
alpha = linspace(1.05,20,10000);

for tMV = 40:10:200

    for s = 10:10:tMV-10
        k = s;
        p = tMV - s;
        
        boundPrecSTE(:,s/10) = sqrt(2) .* (1./sqrt(2.*alpha-1) .* sqrt( (k.^(-2*alpha+1) - n.^(-2*alpha+1)) ) + (k-1)/p .*1./sqrt(alpha-1).*k.^(-alpha/2).*sqrt( (k.^(-alpha+1) - n.^(-alpha+1)) ) + (k.^(-alpha+1) - n.^(-alpha+1)) .* min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p),(k-1)/p) .* 1./(alpha-1));
        boundBIGNys(:,s/10) = ((1 + (k+5-1)/(p)).*1./(alpha-1).*((k+5).^(-alpha+1) - n.^(-alpha+1)));
        badboundPrecSTE(:,s/10) = sqrt(2) .* (1./sqrt(2.*alpha-1).* sqrt( (k.^(-2*alpha+1) - n.^(-2*alpha+1)) ) + k.^(-alpha+1) .* (k-1)/p .* 1./(alpha-1));
        
        % l = 5; %Lanczos its
    
        
        % figure(s/10)
        % semilogy(alpha, boundPrecSTE,'r')
        % hold on
        % semilogy(alpha, badboundPrecSTE,'b')
        % hold on
        % semilogy(alpha, boundBIGNys,'g')
        % hold on
        % semilogy(alpha, boundLanczos,'-m')
        % hold on
        % plot(coef,((1 + (k+3-1)/(p+3)).*1./(coef-1).*(k+3).^(-coef+1) ),'r*')
        % hold on
        % plot(badcoef,((1 + (k+3-1)/(p+3)).*1./(badcoef-1).*(k+3).^(-badcoef+1) ),'b*')
        % legend('bound PrecSTE', 'worst bound PrecSTE', 'bound BIGNystrom', 'bound Lanczos 6')
        % xlabel('$\alpha$', 'Interpreter', 'latex')
        % ylabel('errors')
        % title('Bounds for the estimates for algebraic decay', 'Interpreter','latex')
    
    end
    
    BESTPrecSTE = min(boundPrecSTE');
    BESTBIGNys = min(boundBIGNys');
    BESTbadPrecSTE = min(badboundPrecSTE');
    
    tol = 1e-3;
    indexbound = abs(BESTPrecSTE - BESTBIGNys)./abs(BESTBIGNys) < tol;
    intersection = find(indexbound);
    
        if ~isempty(intersection) 
            coef = alpha(intersection(1));
            nu(tMV/10-3) = coef(1);
        else 
            coef = alpha(end);
            nu(tMV/10-3) = coef(1);
        end
        
        badindexbound = abs(BESTbadPrecSTE - BESTBIGNys)./abs(BESTBIGNys) < tol;
        badintersection = find(badindexbound);
    
        if ~isempty(badintersection) 
            badcoef = alpha(badintersection(1));
            badnu(tMV/10-3) = badcoef(1);
        else 
            badcoef = alpha(end);
            badnu(tMV/10-3) = badcoef(1);
        end
    
    figure(1)
    semilogy(alpha, BESTPrecSTE,'r')
    hold on
    semilogy(alpha, BESTbadPrecSTE,'b')
    hold on
    semilogy(alpha, BESTBIGNys,'g')
    
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
