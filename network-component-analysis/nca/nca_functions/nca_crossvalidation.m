function [dependent_operons,X_pred,X_true,dum1,dum0] = nca_crossvalidation(W,X,time_training)

% [dependent_operons,X_pred,X_true,X_pred_all,X_true_all] = nca_crossvalidation(W,X,time_training)
%
% perform crossvalidation of NCA results
% 1. choose a subset of operons necessary to estimate the TF time courses    
%    the time courses of the remaining genes ('dependent_operons') are predicted from the model

dependent_operons=[];
X_pred = [];
X_true = [];

A_hyp = W.data .* randn(size(W.data));

% find dependent operons
%[echelon,independent_operons] = rref(A_hyp');
%T                    = echelon(1:length(independent_operons),:);
%dependent_operons    = setdiff(1:size(A_hyp,1),independent_operons);

fullrank = rank(A_hyp);
dependent_operons = [];
for it = 1:size(X,1),
  if rank(A_hyp(setdiff(1:size(A_hyp,1),it),:)) == fullrank,  
    dependent_operons=[dependent_operons it]; 
  end
end

% ---

if ~exist('time_training','var'),time_training =  1:20:size(X,2); end

time_test     = setdiff(1:size(X,2),time_training);

[A2,B2] = nca(X(:,time_training),W.data,W.signs);
[A2,B2] = nca_adjust_signs(A2,B2,W);

X_true = X(:,time_test);
X_pred = zeros(size(X,1),length(time_test));

for it = 1:length(dependent_operons),

fprintf('Crossvalidation: predicting operon %d out of %d\n',it,length(dependent_operons));
choose = setdiff(1:size(A_hyp,1),dependent_operons(it));
B1 = pinv(A2(choose,:))*X(choose,time_test);
X_pred(it,:) = A2(dependent_operons(it),:) * B1;
end

if nargout>3,

dum0 = X(:,time_test); dum0(dependent_operons,:) = X_true;
dum1 = zeros(size(X,1),size(X_true,2)); dum1(dependent_operons,:) = X_pred;

end
