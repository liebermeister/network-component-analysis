function [W_new,A_new,changes,X_pred] = modify_signs(X,W,A,B,W_signs,X_error,p_sign_error,fixed_TF,verbose,general,lambda_dep,regularise_flag)

% [W_new,A,B,changes,X_pred] = relax_signs(X,W,B,W_signs,p_sign_error,verbose)
%
% updates the sign restrictions in W based on a matrix of (supposedly) known signs
% and an error probability. 
% Note that the matrix B is kept constant during the optimisation!

lambda_regular_A = 0;

if ~exist('verbose','var'),  verbose = 0; end 
if ~exist('p_sign_error','var'),  p_sign_error = 0.05; end 
if ~exist('W_signs','var'),  W_signs = W.signs; end 
if size(X_error)==[1 1], X_error = X_error * ones(size(X)); end 

n_genes = size(X,1);

signs_new = zeros(size(W.signs));

for gene = 1:n_genes,

 if verbose,  fprintf('Checking gene %d\n',gene);end
 
  x = X(gene,:);
  x_error = X_error(gene,:);
  TF_ind = find(W.data(gene,:));
  TF_ind_restr = TF_ind(W_signs(gene,TF_ind)~=0);
  if isempty(TF_ind_restr), signs_new(gene,:) = W_signs(gene,:);
      else, 
  fxd_sgns = W_signs(gene,TF_ind_restr);
  a = A(gene,TF_ind);
  combinations = num2bit(0:2^length(TF_ind_restr)-1);
  clear log_likelihood log_posterior log_prior
   for it = 1:size(combinations,1),
     this_fxd_sgns = W.signs(gene,TF_ind);
     this_fxd_sgns(W_signs(gene,TF_ind)~=0) = fxd_sgns.*combinations(it,:);
     log_prior(it) =     log(p_sign_error) * sum(1-combinations(it,:))...
	             + log(1-p_sign_error) * sum(combinations(it,:));
     A_new  = optimise_A_given_B(x,ones(size(this_fxd_sgns)),B(TF_ind,:),this_fxd_sgns,lambda_regular_A,a,general,lambda_dep,regularise_flag);
     x_pred = A_new*B(TF_ind,:);
     log_likelihood(it) = nanmean(-(x-x_pred).^2./(2*x_error.^2),2);
     log_posterior(it) = log_prior(it) +  log_likelihood(it);
  end
    if verbose,
      W_signs(gene,TF_ind)
      a
%      combinations *diag(fxd_sgns)
      [      log_prior'      log_likelihood'      log_posterior']
      pause
  end

  [dum,best] = min(-log_posterior);
  signs_new(gene,TF_ind_restr) =  W.signs(gene,TF_ind_restr) .* combinations(best,:);
  end
end

signs_new(:,fixed_TF) = W.signs(:,fixed_TF);

W_new = W;
W_new.signs = signs_new;
A_new = optimise_A_given_B(X,W_new.data,B,W_new.signs,lambda_regular_A,A,general,lambda_dep,regularise_flag);

dummy = sparse(abs(W.signs) - abs(W_new.signs));
[i1,i2] = ind2sub(size(W.data),find(full(dummy)));
changes = [i2 i1];
 
if nargout ==5, X_pred = A_new*B; end
