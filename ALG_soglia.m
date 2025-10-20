% rho = 0.55;
% rho = 0.7;
% rho = 0.85;
rho = 1;

alpha = -rho*2.^linspace(1,9,9);

kalpha = zeros(9,1); 
kappa = zeros(8,1);
m = 3;

for j=1:9

    lhs = (alpha(j)+1)*(2*alpha(j)+1)*(2*(m+1)); %3 = (m+1)/2
    lhs - 0.25*(2*(m+1))*(8*alpha(j).^2+12*alpha(j)+4) 
    (2*alpha(j)+1)*(4*alpha(j)+1)*(2*(m+1)) - (2*(m+1))*(8*alpha(j).^2+6*alpha(j)+1) 
    for k = 5:1:50000
        rhs = ((k+2+m)/(k))^(2*alpha(j)) * (k+m+2)^2;
        if rhs>=lhs
            kalpha(j) = k;
            break
        end
    end
    
    if j>1
        kappa(j-1) = kalpha(j) - 2*kalpha(j-1);
    end

end

