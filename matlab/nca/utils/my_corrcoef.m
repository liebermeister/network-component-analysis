function [c, pvalue_one_tailed, pvalue_two_tailed] = my_corrcoef(A,B)

% [c, pvalue_one_tailed, pvalue_two_tailed] = my_corrcoef(A,B)
%
% A, B row vectors

c = corrcoef([A' B']);
c = c(1:size(A,1), end-size(B,1)+1:end);

if nargout>1,
if size(A,1)==1,
    if abs(c) ~= 1,
n = size(A,2);
t = c * sqrt((n-2)/(1-c^2));
pvalue_one_tailed = 1 - tcdf(t,n-2);
pvalue_two_tailed = 2 * pvalue_one_tailed;
end
end
end