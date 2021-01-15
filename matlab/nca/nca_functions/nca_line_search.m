function [A_new,B_new,lambda] = nca_line_search(X,A,B,A_new,B_new);

delta_A = A_new - A;
delta_B = B_new - B;

lambda = fminbnd(@error,1,2,[],X,A,B,delta_A,delta_B);

A_new = A + lambda * delta_A;
B_new = B + lambda * delta_B;

function e = error(lambda,X,A,B,delta_A,delta_B)

 e = nanmedian(nanmedian( (X - (A+lambda *delta_A) * (B+lambda*delta_B)).^2)');
