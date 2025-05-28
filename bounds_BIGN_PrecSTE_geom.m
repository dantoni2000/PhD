clear, close all

n = 1000;
alpha = linspace(0,0.95,10000); 

for k = 10:10:100

    p = k;
    
    eqboundPrecSTE = sqrt(2) * (sqrt(1-alpha.^(2*n)).*sqrt(1-alpha) + sqrt(1-alpha.^(n)).*(k-1)/p .* sqrt(1-alpha.^2) + min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p), (k-1)/p) .* (1-alpha.^n) .*sqrt(1+alpha));
    
    eqboundBIGNys = alpha.^3 .* (1+(k+3-1)/(p+3)) .*(1-alpha.^n).* sqrt(1+alpha);
    
    eqbadboundPrecSTE = sqrt(2) * (sqrt(1-alpha.^(2*n)).*sqrt(1-alpha) + (k-1)/p * (1-alpha.^n).* sqrt(1+alpha));
    
    normE = alpha.^(k).*(1+2*k/(p-1) + 2*exp(2)*(k+p)/(p^2-1)./(1-alpha)); 
    boundLanczos = 2 .* (sqrt(normE+2) +1) .* (log(2) + normE) .* n .* ((sqrt(normE +2)-1)./((sqrt(normE +2)+1))).^12;

    % l = 4; %Lanczos its
    
    tol = 1e-3;
    indexbound = abs(eqboundPrecSTE - eqboundBIGNys)./abs(eqboundBIGNys) < tol;
    intersection = find(indexbound);

    if ~isempty(intersection) 
        coef = alpha(intersection(1));
        gamma(k/10) = log(1/coef(1));
    else 
        coef = alpha(end);
        gamma(k/10) = log(1/coef(1));
    end
    
    badindexbound = abs(eqbadboundPrecSTE - eqboundBIGNys)./abs(eqboundBIGNys) < tol;
    badintersection = find(badindexbound);

    if ~isempty(badintersection) 
        badcoef = alpha(badintersection(1));
        badgamma(k/10) = log(1/badcoef(1));
    else 
        badcoef = alpha(end);
        badgamma(k/10) = log(1/badcoef(1));
    end
    
    
    figure(k/10)
    plot(alpha, eqboundPrecSTE,'r')
    hold on
    plot(alpha, eqbadboundPrecSTE,'b')
    hold on
    plot(alpha, eqboundBIGNys,'g')
    % hold on
    % plot(alpha, alpha.^4 .* (1+(k+4-1)/(p) .* sqrt(1+alpha)),'k')
    hold on
    plot(coef,coef^3 .* (1+(k+3-1)/(p+3)) .*(1-coef.^n).* sqrt(1+coef),'r*')
    hold on
    plot(badcoef,badcoef^3 .* (1+(k+3-1)/(p+3)) .*(1-badcoef.^n).* sqrt(1+badcoef),'b*')
    legend('reduced bound PrecSTE', 'reduced worst bound PrecSTE', 'reduced bound BIGNystrom')
    xlabel('$\alpha$', 'Interpreter', 'latex')
    ylabel('errors')
    title('Bounds for the estimates for geometric decay', 'Interpreter','latex')
    

    % true bounds
    boundPrecSTE = sqrt(2) * alpha.^k .* (sqrt(1-alpha.^(2*n))./sqrt(1-alpha.^2) + (k-1)/p .*sqrt(1-alpha.^n)./ sqrt(1-alpha) + min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p), (k-1)/p) .*(1-alpha.^n)./ (1-alpha));
    boundBIGNys = alpha.^(k+3) .* (1+(k+3-1)/(p+3)) .*(1-alpha.^n)./(1-alpha);
    badboundPrecSTE = sqrt(2) * alpha.^k .*(sqrt(1-alpha.^(2*n))./sqrt(1-alpha.^2) + (k-1)/p .*(1-alpha.^n)./ (1-alpha));
    
    figure(k/10 + 50)
    semilogy(alpha, boundPrecSTE,'r')
    hold on
    semilogy(alpha, badboundPrecSTE,'b')
    hold on
    semilogy(alpha, boundBIGNys,'g')
    hold on
    semilogy(alpha, boundLanczos,'-m')
    % hold on
    % plot(alpha, alpha.^4 .* (1+(k+4-1)/(p) .* sqrt(1+alpha)),'k')
    hold on
    plot(coef,coef.^(k+3) .* (1+(k+3-1)/(p+3)) .*(1-coef.^n)./(1-coef),'r*')
    hold on
    plot(badcoef,badcoef.^(k+3) .* (1+(k+3-1)/(p+3)) .*(1-badcoef.^n)./(1-badcoef),'b*')
    legend('bound PrecSTE', 'worst bound PrecSTE', 'bound BIGNystrom', 'bound Lanczos 6')
    xlabel('$\alpha$', 'Interpreter', 'latex')
    ylabel('errors')
    title('Bounds for the estimates for geometric decay', 'Interpreter','latex')
end

figure(100)
% plot(10:10:100,gamma,'-*r')
% hold on
% plot(10:10:100,badgamma,'-*b')
% hold on
plot(10:10:100,max(gamma,badgamma),'-*k')
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\gamma$', 'Interpreter', 'latex')
title('Coefficient $\gamma$ for which $trace \log((I + e^{-\gamma x_{ii}}))$ can be better approximated with PrecSTE', 'Interpreter','latex')
