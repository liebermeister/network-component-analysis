function derivative = differentiator(data,M)

% numerical derivative using hamming window
% problems: 
%  1 artefacts at the boundaries,
%  2 number of data points decreases

if ~exist('M','var'), M   = 31;  end

alpha = (M-1)/2;
n=0:M-1;

hd = cos(pi*(n-alpha))./(n-alpha); hd(alpha+1)=0;
w_ham = (hamming(M))';
h = hd .* w_ham;
derivative = conv2(conv2(data,h,'same'),w_ham,'valid');
