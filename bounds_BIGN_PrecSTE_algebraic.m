clear, close all

n = 1000;
alpha = linspace(1.05,20,10000);

for k = 10:10:100

    p = k; % p = 25;
    
    % boundPrecSTE = sqrt(2) .* ((1./sqrt(2.*alpha-1) + (k-1)/p .*1./sqrt(alpha-1)) .* k.^(-alpha+0.5) + k.^(-alpha+1) .* exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p) .* 1./(alpha-1));
    boundPrecSTE = sqrt(2) .* (1./sqrt(2.*alpha-1) .* sqrt( (k.^(-2*alpha+1) - n.^(-2*alpha+1)) ) + (k-1)/p .*1./sqrt(alpha-1).*k.^(-alpha/2).*sqrt( (k.^(-alpha+1) - n.^(-alpha+1)) ) + (k.^(-alpha+1) - n.^(-alpha+1)) .* min(exp(1) .* sqrt(k+p)/p .* sqrt((k-1)/p),(k-1)/p) .* 1./(alpha-1));
    % boundBIGNys = ((1 + (k+3-1)/(p+3)).*1./(alpha-1).*(k+3).^(-alpha+1) );
    boundBIGNys = ((1 + (k+3-1)/(p+3)).*1./(alpha-1).*((k+3).^(-alpha+1) - n.^(-alpha+1)));
    % badboundPrecSTE = sqrt(2) .* ((1./sqrt(2.*alpha-1)) .* k.^(-alpha+0.5) + k.^(-alpha+1) .* (k-1)/p .* 1./(alpha-1));
    badboundPrecSTE = sqrt(2) .* (1./sqrt(2.*alpha-1).* sqrt( (k.^(-2*alpha+1) - n.^(-2*alpha+1)) ) + k.^(-alpha+1) .* (k-1)/p .* 1./(alpha-1));
    
    normE = k.^(-alpha).*(1+2*k/(p-1) + 2*exp(2)*(k+p)/(p^2-1)*k./(alpha-1)); 
    boundLanczos = 2 .* (sqrt(normE+2) +1) .* (log(2) + normE) .* n .* ((sqrt(normE +2)-1)./((sqrt(normE +2)+1))).^12;
    % l = 4; %Lanczos its

    tol = 1e-3;
    indexbound = abs(boundPrecSTE - boundBIGNys)./abs(boundBIGNys) < tol;
    intersection = find(indexbound);

    if ~isempty(intersection) 
        coef = alpha(intersection(1));
        nu(k/10) = coef(1);
    else 
        coef = alpha(end);
        nu(k/10) = coef(1);
    end
    
    badindexbound = abs(badboundPrecSTE - boundBIGNys)./abs(boundBIGNys) < tol;
    badintersection = find(badindexbound);

    if ~isempty(badintersection) 
        badcoef = alpha(badintersection(1));
        badnu(k/10) = badcoef(1);
    else 
        badcoef = alpha(end);
        badnu(k/10) = badcoef(1);
    end

    
    figure(k/10)
    semilogy(alpha, boundPrecSTE,'r')
    hold on
    semilogy(alpha, badboundPrecSTE,'b')
    hold on
    semilogy(alpha, boundBIGNys,'g')
    hold on
    semilogy(alpha, boundLanczos,'-m')
    hold on
    plot(coef,((1 + (k+3-1)/(p+3)).*1./(coef-1).*(k+3).^(-coef+1) ),'r*')
    hold on
    plot(badcoef,((1 + (k+3-1)/(p+3)).*1./(badcoef-1).*(k+3).^(-badcoef+1) ),'b*')
    legend('bound PrecSTE', 'worst bound PrecSTE', 'bound BIGNystrom', 'bound Lanczos 6')
    xlabel('$\alpha$', 'Interpreter', 'latex')
    ylabel('errors')
    title('Bounds for the estimates for algebraic decay', 'Interpreter','latex')

end

figure(100)
% plot(10:10:100,nu,'-*r')
% hold on
% plot(10:10:100,badnu,'-*b')
% hold on
plot(10:10:100,max(nu,badnu),'-*k')
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\nu$', 'Interpreter', 'latex')
title('Coefficient $\nu$ for which $trace \log((I + x_{ii}^{-\nu}))$ can be better approximated with PrecSTE', 'Interpreter','latex')
