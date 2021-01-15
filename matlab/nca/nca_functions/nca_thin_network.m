function [W_red,A_red,n_omit,n_keep,ssr_change] = nca_thin_network(X,A,B,X_pred,percentage,W)

% DIRECT increase in sum of squared residuals if a single influence value is set to 0 
%
% for influence Aij:
%    || X_pred_i - (X_i - Aij * Bj) ||^2 - || X_pred_i - X_i ||^2
% =  || (X_pred_i - X_i) - Aij * Bj ||^2 - || X_pred_i - X_i ||^2
% =  -2 < (X_pred_i - X_i),  Aij * Bj > + < Aij * Bj, Aij * Bj >
%
% yielding A.^2 * diag(diag(B*B')) -  2 * A .* ( (X_pred-X) * B') 
%
% percentage: approximate percentage by which the SSR may increase

indfin = find(isfinite(X));
ssr_tot    = sum(sum((X_pred(indfin)-X(indfin)).^2));
ssr_change =  A.^2 * diag(diag(B*B')) -  2 * A .* ( (X_pred-X) * B');

ssr_change(ssr_change<0)==0; 
% actually, ssr_change should never be negative, but slightly 
% negative values occur ... (incomplete nca optimisation?)

sorted_ssr_changes = sort(ssr_change(find(W.data)));
n_omit = sum(cumsum(sorted_ssr_changes) < percentage*ssr_tot);
threshold = sorted_ssr_changes(n_omit);

A_red = A;
A_red(ssr_change < threshold) = 0;
W_red = W;
W_red.data(ssr_change < threshold) = 0;
W_red.signs = sign(A_red);
n_keep = sum(W.data(:))-n_omit;
