function [ truncval ] = truncNorm(x,mu,vars,a,b)
% this will input a point where you would like to know the value of a
% truncated normal pdf with the truncation points of a and b

if nargin < 1, x = [0.5 0.75]; end
if nargin < 2, mu = 2.4; end
if nargin < 3, vars = 0.8; end
if nargin < 4, a = 0; end
if nargin < 5, b = Inf; end

if iscell(x), x = cell2mat(x); end
stdx = sqrt(vars);
truncval = normpdf(x,mu,stdx)./(normcdf((b-mu)./stdx)-normcdf((a-mu)./stdx));
idx = x < a;
truncval(idx) = 0;

end