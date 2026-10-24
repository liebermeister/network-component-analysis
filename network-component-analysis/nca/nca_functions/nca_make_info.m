function info = nca_make_info(A_list,B_list,Wdata,Wsigns,X,lambda_A,force_B1,general,lambda_dep,regularise_flag)

if ~exist('lambda_A','var'), lambda_A = 0; end
if ~exist('force_B1','var'), force_B1 = []; end

info = [];

A_sum=zeros(size(A_list{1}));
B_sum=zeros(size(B_list{1}));
A_sqr=zeros(size(A_list{1}));
B_sqr=zeros(size(B_list{1}));

for it = 1:length(A_list),
  sgns = diag(sign(B_list{it}*B_list{1}'));
  sgns(sgns==0) = 1;
  sgns(find((sgns==-1)' .* (sum(abs(Wsigns))))) = 1;
  A_list{it} = A_list{it}*diag(sgns);
  B_list{it} = diag(sgns)*B_list{it};
  
  A_sum = A_sum+A_list{it};
  A_sqr = A_sqr+A_list{it}.^2;
  B_sum = B_sum+B_list{it};
  B_sqr = B_sqr+B_list{it}.^2;
  errors(it) = nanmean(nanmean( (X - A_list{it} * B_list{it}).^2)');
end

A.mean = A_sum/length(A_list);
B.mean = B_sum/length(A_list);
A.std = sqrt(A_sqr/length(A_list)-A.mean.^2);
B.std = sqrt(B_sqr/length(B_list)-B.mean.^2);
A.list = A_list;
B.list = B_list;
[dum,best] = min(errors);
A.best = A.list{best};
B.best = B.list{best};

weights =  exp( - (errors-min(errors)) / ( 0.5*(median(errors) - min(errors))+10^-10));

dum=[]; for it=1:length(B.list),
  dum = [dum; weights(it)*reshape(B.list{it},1,prod(size(B.mean)))]; 
end; 

B.robust = reshape(sum(dum,1)/sum(weights),size(B.mean,1),size(B.mean,2));
A.robust = optimise_A_given_B(X,Wdata,B.robust,Wsigns,lambda_A,A.best,general,lambda_dep,regularise_flag);
[A.robust,B.robust] = nca_rescale(A.robust,B.robust,force_B1);

X_pred = A.robust*B.robust;
X_pred(isnan(X)) = nan;

info.A = A;
info.B = B;
info.errors = errors;
