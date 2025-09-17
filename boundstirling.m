clear, close all

for x = 1:1:20
    for m = 2:1:100
    
        a = (m-1)/2;
        b = 3/2;
        y = 2e-6*2^x;
        % y = 2^x;
        z = y/(2*m);
    
        f = @(u)(exp(-z.*(u./(1-u))).*u.^(a-1).*(1-u).^(b-1));
        true(m-1,x) = 1/beta(a,b) * integral(f,0,1);
        gg(m-1,x) = sqrt(m)/sqrt(m+1) * (1+sqrt(y))*exp(-sqrt(y));
        
        value = sqrt((a)/(z));
        
        ggg(m-1,x) = exp(-z*value)*(betainc(value/(1+value),a,b));
        g(m-1,x) = exp(-z*2*a);
        s=@(x)(sqrt(m./(m+1))*1./(2*sqrt(x)).*exp(-x.*(1-1./(2*m))));
        r=@(x)(sqrt(m)/sqrt(m+1) *(1+sqrt(x))./(2*sqrt(x)).*exp(-x/2-sqrt(x)));
        integral(s,0,2)
        sqrt(m*pi)/(2*sqrt((m+1)*(1-1/(2*m))))*erf(sqrt((1-1/(2*m))*2))
        integral(r,0,2)

  %       bah = @(u) ( ...
  %   (u.^(m/2 - 3/2).*exp((u.*x)./(2*m.*(u - 1))).*log(u).*(1 - u).^(1/2)) ...
  %       /(2*beta(3/2, m/2 - 1/2)) ...
  % + (u.^(m/2 - 3/2).*exp((u.*x)./(2*m.*(u - 1))).*(psi(m/2 + 1)/2 - psi(m/2 - 1/2)/2).*(1 - u).^(1/2)) ...
  %       /beta(3/2, m/2 - 1/2) ...
  % + (u.*u.^(m/2 - 3/2).*x.*exp((u.*x)./(2*m.*(u - 1)))) ...
  %       /(2*m.^2.*beta(3/2, m/2 - 1/2).*(1 - u).^(1/2)) ...
  %   );
  %       dTRUE_int(m-1,x) = integral(bah,0,1);
  %       dTRUE_int(dTRUE_int>0)
    end
end

for j=1:1:20
    figure(j)
    semilogy(true(1:end,j),'g')
    hold on
    semilogy(gg(1:end,j),'b')
    hold on
    semilogy(g(1:end,j),'r')
    hold on
    semilogy(ggg(1:end,j),'k')
end

% for s=1:1:30
%     figure(20 + s)
%     semilogy(true(s,1:end))
%     hold on
%     semilogy(gg(s,1:end))
% end