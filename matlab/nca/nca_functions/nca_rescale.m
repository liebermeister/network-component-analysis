function [A, B] = nca_rescale(A, B, force_B1)

% [A,B] = nca_rescale(A,B,force_B1)
%
% rescale A and B such that the mean square value of B == 1

n_force            = size(force_B1,1);
scales             = sqrt( mean(B(n_force+1:end,:)'.^2)' ) + 0.00000001;
B(n_force+1:end,:) = diag(1./scales) * B(n_force+1:end,:);
A(:,n_force+1:end) = A(:,n_force+1:end) * diag(scales);
