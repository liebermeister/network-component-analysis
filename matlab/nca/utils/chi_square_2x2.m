% function [resultat,pval] = chi_square_2x2(matrix,alpha)
%
% matrix: 2x2-kontingenztafel, nullhyp: variablen stat. unabh.
% resultat: nullhyp verwerfen auf alpha-konfidenzniveau?

function [resultat,pval] = chi_square_2x2(g,alpha)

a=g(1,1);
b=g(1,2);
c=g(2,1);
d=g(2,2);

n   = a+b+c+d;
chi_sq = n*(a*d-b*c)^2/((a+b)*(a+c)*(b+d)*(c+d));

if a+b<6 | c+d < d
  fprintf('Sample too small. Chi^2 test should not be applied\n');
end
%chi_sq_star = (n-1)*(a*d-b*c)^2/((a+b)*(a+c)*(b+d)*(c+d))  % better for small samples

%plot(0:0.1:10,2*(1-CDFnormal(sqrt(0:0.1:10))))
%[0:0.1:10;2*(1-CDFnormal(sqrt(0:0.1:10)))]

pval = 2*(1-CDFnormal(sqrt(chi_sq)));

resultat = pval < alpha;

%---------------------------------------------

function Iu = CDFnormal(u)
t  = 1./(1 + 0.2316419 * u);
Iu = 1 - 1/sqrt(2 * pi) * exp(-1/2 .* u.^2) .* ...
     ( 0.319381530 * t   - 0.356563782 * t .^2 + ...
       1.781477937 * t .^3 - 1.821255978 * t .^4 + ...
       1.330274429 * t .^5) ;

