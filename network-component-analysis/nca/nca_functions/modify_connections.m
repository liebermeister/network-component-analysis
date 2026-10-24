function W_new = modify_connections(X,W,A,B,P,known_signs,sigma_0,p_sign_error,n_threshold_possible,verbose,general,lambda_dep,regularise_flag)

% W_new = choose_new_W(X,W,B,P,known_signs,sigma_0,p_sign_error,n_threshold_possible,verbose)
%
% n_threshold_possible: maximal number of inputs per gene up to which all combinations are checked
% sigma_0: estimated std dev of measurements
% P prior probabilities of for non-zero element of W.data
% p_sign_error: prior probability of a wrong sign in W.signs

if ~exist('verbose','var'),  verbose = 1; end 

lambda_regular_A = 0;
n_genes = size(X,1);

for gene = 1:n_genes,
 
  if verbose,  fprintf('Checking gene %d\n',gene);end
 
  x = X(gene,:);
  prior = P(gene,:);
  possible_inputs = find(prior);
  n_possible = sum(prior~=0);

  logprior_yes = log(P(gene,possible_inputs));
  logprior_no  = log(1-P(gene,possible_inputs));
 
  if n_possible <= n_threshold_possible,
    combinations = num2bit(0:2^n_possible-1);
  else 
    current = (W.data(gene,:)~=0);
    current = current(possible_inputs);  
    combinations = [current; mod(repmat(current,length(current),1)+eye(length(current)),2)];
  end
  clear log_prior log_likelihood log_posterior TF_choice

  for it = 1:size(combinations,1),
       log_prior(it) = sum(logprior_yes(find(combinations(it,:)))) ...
                 + sum(logprior_no(find(1-combinations(it,:))));
       TF_choice{it} = possible_inputs(find(combinations(it,:)));
       if length(TF_choice{it}),
%        B_explain = B(possible_inputs(find(combinations(it,:))),:);
%        errors    = x' - B_explain'*inv (B_explain*B_explain')*B_explain * x';
         A_new  = optimise_A_given_B(x,ones(1,length(TF_choice{it})),B(TF_choice{it},:),...
                                     W.signs(gene,TF_choice{it}),lambda_regular_A,A,...
                                     general,lambda_dep,regularise_flag);
         x_pred = A_new*B(TF_choice{it},:);
       else,
         x_pred = zeros(size(x));
       end
         errors    = x - x_pred;
%      log_likelihood(it) = sum(-(x-x_pred).^2./(2*x_error.^2));
         if nanmean(errors.^2,2) > 10^20, 
           A_new
         end
       log_likelihood(it) = -1/(2*sigma_0^2)*nanmean(errors.^2,2);
       log_posterior(it) = log_prior(it) + log_likelihood(it);
  end

  if verbose, 
    [dum,order] = sort(-log_posterior);
    for it = 1:length(order),
      fprintf('%s: %f %f %f\n',num2str(TF_choice{order(it)}), log_prior(order(it)),...
              log_likelihood(order(it)),log_posterior(order(it)));
    end
  end

  [dum,best] = min(-log_posterior);
  best_choice{gene} = TF_choice{best};

  if verbose, fprintf('Best choice: '); fprintf('%d ', best_choice{gene}); fprintf('\n'); end

end

W_new = W;

for it = 1:length(best_choice), W_new.data(it,best_choice{it}) = 1; end

unused_TF = find(1-(sum(W_new.data,1)>0));

for it =1:length(unused_TF),
 [dum,gg] = max(P(:,unused_TF(it)));
 W_new.data(gg,unused_TF(it))=1;
end

