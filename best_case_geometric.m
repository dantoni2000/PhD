clear, close all

n = 1000;
alpha = linspace(0.69,0.94,5000); 

for tMV = 10:10:90

    for s = 4:2:tMV-4
        k = s;
        p = tMV - s;
        
        boundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k) .* (sqrt(1-alpha.^(2*n))./sqrt(1-alpha.^2) + (k-1)/p .*sqrt(1-alpha.^n)./ sqrt(1-alpha) + min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p), (k-1)/p) .*(1-alpha.^n)./ (1-alpha));
        boundBIGNys(:,s/2-1) = alpha.^(k+5) .* (1+(k+5-1)/(p)) .*(1-alpha.^n)./(1-alpha);
        badboundPrecSTE(:,s/2-1) = sqrt(2) * alpha.^(k) .*(sqrt(1-alpha.^(2*n))./sqrt(1-alpha.^2) + (k-1)/p .*(1-alpha.^n)./ (1-alpha));
        
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
