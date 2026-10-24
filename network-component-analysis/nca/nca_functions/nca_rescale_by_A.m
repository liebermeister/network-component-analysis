function [A,B] = nca_rescale_by_A(A,B,Wdata);

% scaling such that the mean absolute value of relevant A elements is 1

dd = sum(abs(A))./sum(Wdata);
 A = A * diag(1./dd);
 B = diag(dd) * B;