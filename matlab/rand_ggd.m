function x = rand_ggd(n,N,mu,beta,rho)
% Generate N Generalized Gaussian m-dimensional vectors with means mu(:),
% inverse scales beta(:), and parameters rho(:).
% EXAMPLE:
%   N = 10e3;
%   mu = 0;
%   alpha = 1.27;
%   sigma = 1;
%   beta = sigma*sqrt(gamma(1/alpha)/gamma(3/alpha));
%   X = rand_ggd(1,N,mu,beta,alpha);
%
%   Example generates 10e3 GGD samples with shape = 1.27.
%   Script modified from original:
%
%  [1] Palmer, J. (2011, September 14). Generalized Gaussian Probability
%  Density Function. Retrieved from https://sccn.ucsd.edu/wiki/Generalized_Gaussian_Probability_Density_Function
if length(mu) == 1
    mu = mu*ones(1,n);
    beta = beta*ones(1,n);
    rho = rho*ones(1,n);
end
x = zeros(n,N);
for i = 1:n
    x(i,:) = mu(i) + (1/sqrt(beta(i))) * ( mygamrnd(1/rho(i),1,1,N)).^(1/rho(i) ) .* ((rand(1,N)<0.5)*2-1) ;
end

end

function r = mygamrnd(a,b,rows,columns)
% This is essentially copied from the gamrnd function in the Matlab Stats
% Toolbox. If you have that toolbox, you can just use the code above and
% replace mygamrnd() with gamrnd().
% Initialize r to zero.
lth = rows*columns;
r = zeros(lth,1);
a = a(:); b = b(:);
scalara = (length(a) == 1);
if scalara
    a = a*ones(lth,1);
end
scalarb = (length(b) == 1);
if scalarb
    b = b*ones(lth,1);
end
% If a == 1, then gamma is exponential. (Devroye, page 405).
k = find(a == 1);
if any(k)
    r(k) = -b(k) .* log(rand(size(k)));
end
k = find(a < 1 & a > 0);
% (Devroye, page 418 Johnk's generator)
if any(k)
    c = zeros(lth,1);
    d = zeros(lth,1);
    c(k) = 1 ./ a(k);
    d(k) = 1 ./ (1 - a(k));
    accept = k;
    while ~isempty(accept)
        u = rand(size(accept));
        v = rand(size(accept));
        x = u .^ c(accept);
        y = v .^ d(accept);
        k1 = find((x + y) <= 1);
        if ~isempty(k1)
            e = -log(rand(size(k1)));
            r(accept(k1)) = e .* x(k1) ./ (x(k1) + y(k1));
            accept(k1) = [];
        end
    end
    r(k) = r(k) .* b(k);
end
% Use a rejection method for a > 1.
k = find(a > 1);
% (Devroye, page 410 Best's algorithm)
bb = zeros(size(a));
c  = bb;
if any(k)
    bb(k) = a(k) - 1;
    c(k) = 3 * a(k) - 3/4;
    accept = k;
    count = 1;
    while ~isempty(accept)
        m = length(accept);
        u = rand(m,1);
        v = rand(m,1);
        w = u .* (1 - u);
        y = sqrt(c(accept) ./ w) .* (u - 0.5);
        x = bb(accept) + y;
        k1 = find(x >= 0);
        if ~isempty(k1)
            z = 64 * (w .^ 3) .* (v .^ 2);
            k2 = (z(k1) <= (1 - 2 * (y(k1) .^2) ./ x(k1)));
            k3 = k1(find(k2));
            r(accept(k3)) = x(k3);
            k4 = k1(find(~k2));
            k5 = k4(find(log(z(k4)) <= (2*(bb(accept(k4)).*log(x(k4)./bb(accept(k4)))-y(k4)))));
            r(accept(k5)) = x(k5);
            omit = [k3; k5];
            accept(omit) = [];
        end
    end
    r(k) = r(k) .* b(k);
end
r = reshape(r,rows,columns);
end