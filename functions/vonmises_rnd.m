function theta = vonmises_rnd(mu,kappa,sz)
% vonmises_rnd  random angles from a von Mises distribution
%   theta = vonmises_rnd(mu,kappa,sz) returns an array of size sz (as in rand)
%   of angles (rad) drawn from the von Mises distribution of mean mu and
%   concentration kappa (kappa >= 0), with density
%       p(theta) = exp(kappa*cos(theta-mu)) / (2*pi*besseli(0,kappa))
%   Angles are returned in [mu-pi, mu+pi]. Uses the rejection algorithm of
%   Best & Fisher (1979, Applied Statistics 28:152-157) with rand only, so the
%   draws depend on the global random number stream (set with rng).

if nargin < 3
    sz = [1 1];
end
if isscalar(sz)
    sz = [sz sz];
end
if kappa < 0
    error('kappa must be non-negative!')
end

n = prod(sz);
if kappa < 1e-8 % uniform distribution
    theta = mu + pi*(2*rand(sz)-1);
    return
end

% constants of the envelope (wrapped Cauchy)
tau = 1+sqrt(1+4*kappa^2);
rho = (tau-sqrt(2*tau))/(2*kappa);
r   = (1+rho^2)/(2*rho);

% rejection sampling, by batches until n angles are accepted
f = nan(n,1);
nacc = 0;
while nacc < n
    nrem = n-nacc;
    u1 = rand(nrem,1);
    u2 = rand(nrem,1);
    z  = cos(pi*u1);
    fc = (1+r*z)./(r+z);
    c  = kappa*(r-fc);
    accept = (c.*(2-c)-u2 > 0) | (log(c./u2)+1-c >= 0);
    nnew = sum(accept);
    f(nacc+(1:nnew)) = fc(accept);
    nacc = nacc+nnew;
end

% random sign around the mean
u3 = rand(n,1);
theta = mu + (2*(u3 > 0.5)-1).*acos(f);
theta = reshape(theta,sz);

end
