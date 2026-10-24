function [X_shifted,shifts,A_new,X_error_shifted] = modify_shifts(W,X,A,B,sigma_prior,sigma_loglike,X_error,experiments)

% [X_shifted,shifts] = modify_shifts(W,X,B,sigma_prior,sigma_loglike)

% the log likelihood for each gene is computed as: - mean sq. error / (2 * sigma_loglike^2)
% the prior  for each gene is computed as: sq. shift / (2 * sigma_prior^2)

if ~exist('experiments','var'),experiments=[]; end

if length(experiments)>1,
  X_shifted = nan*X;
  shifts = {};
  X_error_shifted = [];
  for it = 1:length(experiments),
    ind = experiments{it};
    ind2 = find(sum(isnan(X(:,ind)),2)==0);
    this_W.data = W.data(ind2,:);
    this_W.signs = W.signs(ind2,:);
    this_B = B(:,ind); this_B(find(isnan(this_B)))=0;


    [this_X_shifted,this_shifts,this_A_new,this_X_error_shifted] = ...
      modify_shifts(this_W,X(ind2,ind),this_B,sigma_prior,sigma_loglike,X_error(ind2,ind));
    X_shifted(ind2,ind) = this_X_shifted;
    shifts = union(shifts,this_shifts);
    X_error_shifted(ind2,ind) = this_X_error_shifted;
  end
  A_new = optimise_A_given_B(X_shifted,W.data,B,W.signs,0,A);

else,

Wdata = W.data;
Wsigns = W.signs;
lambda_regular_A = 0;

B(isnan(B))=0;

if ~exist('sigma_loglike','var'),  sigma_loglike = []; end

if isempty(sigma_loglike),
 A_new  = optimise_A_given_B(X,Wdata,B,Wsigns,lambda_regular_A,A);
 X_pred = A_new * B;
 sigma_loglike =  sqrt(mean(mean((X-A_new*B).^2))); 
end

if ~exist('sigma_prior'), sigma_prior = 5; end

% apply modify shifts to increase quality of fit...

errors = [];

for shift_value = -5:1:5,

 if shift_value <0,
  thisB = B(:,1:end+shift_value+1);
  thisX = X(:,-shift_value:end);
 end
 if shift_value ==0,
  thisB = B;
  thisX = X;
 end
 if shift_value >0,
  thisX = X(:,1:end-shift_value+1);
  thisB = B(:,shift_value:end);
 end

 A_new  = optimise_A_given_B(thisX,Wdata,thisB,Wsigns,lambda_regular_A,A);
 X_pred = A_new * B;
 errors = [errors, mean((thisX-A_new*thisB).^2,2)];

end

errors = errors/(2*sigma_loglike^2);
errors = errors + repmat((-5:1:5).^2/(2*sigma_prior^2),size(errors,1),1);

shift_fine = -5:1:5;
errors_fine = interp1(-5:1:5,errors',shift_fine,'cubic')';

%im(normleft_max(errors_fine-repmat(min(errors_fine')',1,size(errors_fine,2))))

[dum,pos] = min(errors_fine'); 
shifts = shift_fine(pos);

X_shifted = X;
for it =1:size(X,1),
  X_shifted(it,:) = shift_matrix(X(it,:),shifts(it));
  errors_fine(it,:) = shift_matrix(errors_fine(it,:),-shifts(it));
end
%im(normleft_max(errors_fine-repmat(min(errors_fine')',1,size(errors_fine,2))))

X_shifted = fill_nan(X_shifted);
if exist('X_error','var'), 
for it =1:size(X,1),
  X_error_shifted(it,:) = shift_matrix(X_error(it,:),shifts(it));
end

end

A_new  = optimise_A_given_B(X_shifted,Wdata,B,Wsigns,lambda_regular_A,A);

end
