function [A,B] = nca_rescale_by_A(A,B,Wdata);

% scaling -> compromise between std deviation of (B) close to 1
% and mean absolute (relevant) entry of A close to 1

% compute two (competing) requested scalings for A

a = 1./column(sum(abs(A)+10^-10)./sum(Wdata+10^-10));
b = column(std(B'));

% b should be rescaled if it is already larger than 1
%
% to find a good compromise, rescale the vector [a; b] such that 
% it comes as close as possible to [1;1]
% and use the necessary scaling factor for rescaling

lambda = [a+b]./[a.^2+b.^2];

lambda(find([b<1].*[a>b])) = 1./a(find([b<1].*[a>b]));

lambda(lambda==0) = 10^-10;
A = A * diag(1./lambda);
B = diag(lambda) * B;