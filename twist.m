%% Author Information 
% This function was written by Alicia Sit and Francesco Di Colandrea 
% of the University of Ottawa, 2024.
% See http://arxiv.org/abs/2403.00948 for more details.

%% Compute the twist distribution given the tilt distribution
% See Eq.3 of http://arxiv.org/abs/2403.00948.
%
% INPUTS
% beta   : integration parameter
% thetaM : max tilt angle
% theta  : tilt distribution
%
% OUTPUTS
% phi    : twist distribution
function phi = twist(beta, thetaM, theta)
    epPar = 12.0;
    epPerp = 4.0;
    k11 = 6.7E-12;
    k22 = 3.4E-12;
    k33 = 10.6E-12;
    
    gamma = (epPar-epPerp)/epPerp;
    alpha = (k33-k22)/k22;
    kappa = (k33-k11)/k11;
    
    % Eq.3
    fun = @(x) real(beta.*sqrt(1+kappa.*sin(x).^2)./((gFun(x,beta,thetaM).*cos(x).^2).*(1+alpha.*sin(x).^2)) );
    
    phi = zeros(1,length(theta));
    for i=1:length(theta)
        phi(i) = integral(fun,0,theta(i),'arrayvalued', true);
    end
    
    % g function
    function g = gFun(theta,beta,thetaM)
        p1 = (sin(thetaM).^2 - sin(theta).^2)/((1+gamma*sin(theta).^2).*(1+gamma*sin(thetaM).^2));
        p2 = (beta.^2*(1+kappa)/(1+alpha)).*(((1+alpha*sin(thetaM).^2).*cos(thetaM).^2).^(-1) - ((1+alpha*sin(theta).^2).*cos(theta).^2).^(-1));
        g = sqrt(p1+p2);
    end

end