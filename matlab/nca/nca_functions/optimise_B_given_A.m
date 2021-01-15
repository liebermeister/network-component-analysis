function B_new = optimise_B_given_A(X,Wdata,A,lambda_A,lambda_B,force_B1,no_nans,blocks,B,general,lambda_dep,regularise_flag)

% B_new = optimise_B_given_A(X,Wdata,A,lambda_B,force_B1,no_nans,blocks)
%
% X is allowed to contain nans, in a controlled manner: 
%
% presence of nans must be indicated by no_nans = 0, and the data must come in blocks (time series)
% blocks{i} contains the indices of the time points belonging to the ith block
% within each block holds: a gene that contains at least one nan value is considered
% unmeasured.

%if ~exist('no_nans','var'),    no_nans = sum(sum(isnan(X))) ==0; end 
%if ~exist('general','var'),    general = [];                     end 
%if ~exist('lambda_dep','var'), lambda_dep  = 1;                  end 

eval(default('B','[]','force_B1','[]','no_nans','0','blocks','{1:size(X,2)}','general','[]','lambda_dep','0','lambda_A','0','lambda_B','0','regularise_flag','0'));

A_original = A;
X_original = X;

if ~isempty(B),
  goal_old   = nca_goal(X,A,B,general,regularise_flag,lambda_dep,lambda_A,lambda_B);
end

%B_undetermined = (Wdata'*isfinite(X)==0);

% force_B1 have to be the first rows of B

if ~isempty(force_B1),
  n_force = size(force_B1,1);
  X       = X - A(:,1:n_force) * force_B1;
  A       = A(:,n_force+1:end);
  general = setdiff(general,1:n_force)-n_force;
end

if regularise_flag,
  dumX  = zeros(size(X));
  if length(lambda_dep)  == 1,
    dumA = lambda_dep * A; 
    dumA(:,general) = 0;
  else,
    dumA  =  A * diag(lambda_dep); 
  end
  X      = [X; dumX];
  A      = [A; dumA];
end

if no_nans,
  if sum(isfinite(A(:))==0), warning('A matrix is broken!!'),  A(find(~isfinite(A))) = 0; end
  B_new = pinv( A' * A + lambda_B ) * A' * X;
else,
  B_new = [];
  for it = 1:length(blocks),
% ind_vg: indices of genes without nans
    ind_vg = find(sum(isnan(X(:,blocks{it})),2)==0);
    B_new = [B_new, pinv( A(ind_vg,:)' * A(ind_vg,:) + lambda_B ) * A(ind_vg,:)' * X(ind_vg,blocks{it})];
  end
end

if ~isempty(force_B1), B_new = [force_B1; B_new]; end

%B_new        = B_new .* (1-B_undetermined);

goal_new     = nca_goal(X_original,A_original,B_new,general,regularise_flag,lambda_dep,lambda_A,lambda_B);

if ~isempty(B),
if goal_new>goal_old, fprintf('Warning, B matrix gets worse!\n');  end
end