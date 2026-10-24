function [A,B,X_pred,info,A_list,B_list] = nca(...
    X,W_data,W_signs, s1, v1, s2, v2, s3, v3, ...
    s4, v4, s5, v5, s6, v6, s7, v7, s8, v8, ...
    s9, v9, s10, v10, s11, v11, s12, v12, ...
    s13, v13, s14, v14, s15, v15,s16, v16,...
    s17, v17, s18, v18, s19, v19, s20, v20, s21, v21,...
    s22, v22, s23, v23, s24, v24, s25, v25, s26, v26);

%[A, B, X_pred, info, A_list, B_list] = nca(X, W_data, W_signs, more parameters...)
%
%Network component analysis
%
%NCA model: X = A * B + eta, where eta is a matrix of independent standard Gaussian noise
%Determine matrices A (edge weights) and B (TF activities) by maximum likelihood estimation
%
%For creating the data matrix X from experimental data, type 'help nca_expression'  
%For creating the matrices W_data and W_signs, type 'help nca_networks' 
%
%
%USAGE EXAMPLE
% [A, B, changes] = nca(X,W.data,W.signs,'verbose',0);
%
%INPUTS
% X       data matrix
% W_data  matrix (elements    0,1) indicating the positions of non-zero elements of A
% W_signs matrix (elements -1,0,1) indicating known signs of elements of A. unknown sign -> 0
%
%This is a wrapper function that sets the parameters. 
%The actual calculation is done by the function nca_loop.m
%
%In the function call, parameter arguments must come in pairs, parameter name + value
%The following parameters can be used
% 
%ITERATIONS
% initial_choice  'random'  methods to initialise matrix A. 
%                           'random' (default), 'svd', 'given' (requires parameters A_first, B_first)
% repeat           1        number of repetitions of whole estimation (best fit is used in the end)
% repeat_noise_std 0        level of random noise to be added to data in the repetitions
% n_itmax          500      maximal number of iterations
% lambda           1        relaxation factor for the iterative optimisation 
% epsilon          0.0001   convergence threshold for changes of matrix B
% n_accept         repeat   number of (best) repetitions to be reported
%
%MODEL STRUCTURE AND REGULARISATION
% general          [];      indices of general components activating all genes
% A_first                   initial guess for A
% B_first                   initial guess for B
% force_B1         []       first rows of B to be kept fixed
% lambda_dep       0        regularisation coefficient for resolving ambiguities
%                           0 -> no regularisation; scalar -> regularisation for general components
%                           column vector (length #TF) -> individual value for each transcription factor
% lambda_A         0        regulatisation factor for A
% lambda_B         0        regulatisation factor for B
% flag_X_positive  1        omit negative values in the predicted X
%                  
%MISSING VALUES    
% no_nans          1        flag: does nca have to account for missing values in the data?
% nan_blocks       {}       .. if so: list of time points sets (index vectors)  
%                  
%PROGRESS REPORT   
% verbose          1        flag
% graphics_flag    0        flag
% TF_names         {}       list of transcription factor names (for graphics only)


% ----------------------
% default parameter values

n_itmax            = 500; % maximal number of iterations
lambda             = 1;   % relaxation factor for the iterative optimisation 
lambda_A           = 0.;  % regulatisation factor for A
lambda_B           = 0;   % regulatisation factor for B
A_first            = [];  % initial guess for A
B_first            = [];  % initial guess for B
graphics_flag      = 0;
verbose            = 1;
force_B1           = [];
initial_choice     = 'random';
repeat             = 1;
n_accept           = nan;
repeat_noise_std   = 0;
epsilon            = 0.0001; % convergence threshold for B
no_nans            = sum(sum(isnan(X)))==0;
nan_blocks         = {1:size(X,2)};
TF_names           = cellstr(repmat(' ',size(W_data,1),1));
general            = [];
lambda_dep         = 0;
flag_X_positive    = 0; % omit negative values in the predicted X
no_nan_predictions = 1;
n_processes        = 1;
dirname            = [];


% ----------------- input parameters

if(rem(nargin-3,2)==1);  error('Optional parameters should always go by pairs'); end
for i=1:(nargin-3)/2,
  str_param = eval (['s' int2str(i)]);
  val_param = eval (['v' int2str(i)]);
  switch str_param,   
    case {'graphics_flag','verbose','force_B1','initial_choice','n_itmax','lambda',...
          'lambda_A','lambda_B','A_first','B_first','repeat','n_accept','repeat_noise_std','epsilon','no_nans',...
          'nan_blocks','TF_names','general','lambda_dep','no_nan_predictions',...
          'n_processes','dirname'},         
      eval([str_param ' = val_param;']);
  end
end
if isnan(n_accept), n_accept = repeat; end

% -----------------

if ~exist('W_signs','var'), W_signs = zeros(size(W_data));  end 

regularise_flag = 0;
if length(lambda_dep)>1, 
  regularise_flag = 1;
  if size(lambda_dep,2)>1, lambda_dep = lambda_dep'; end
else, 
  if lambda_dep ~= 0, regularise_flag = 1; end
end

% ----------------

if verbose,
  fprintf('Starting NCA\n');
  if length(lambda_dep)>1, 
    display(sprintf(' Individual regularisation terms for TF contributions'));
  else,
    switch lambda_dep,
      case 0,    display(sprintf(' No regularisation terms for TF contribution'));
      otherwise, display(sprintf(' Regularisation term %f for specific TF contributions',lambda_dep));
    end
  end
  switch lambda_A,
    case 0,      fprintf(' No direct regularisation on A\n');
    otherwise,   fprintf(' Regularisation on A with lambda=%f\n', lambda_A);
  end
  switch lambda_B,
    case 0,      fprintf(' No direct regularisation on B\n');
    otherwise,   fprintf(' Regularisation on B with lambda=%f\n', lambda_B);
  end
  switch initial_choice,
    case 'svd',       fprintf(' Initialising the matrices by singular value decomposition\n');
    case 'random',    fprintf(' Initialising the matrices by random guess\n');
    case 'given',     fprintf(' Initialising the matrices by given matrices A and B\n');
  end
end

% -------------------------------

if min(sum(W_data')) == 0 | min(sum(W_data)) == 0,   % if there are irrelevant genes or TF 
  
  warning('Warning: irrelevant genes or TF. Results may be incorrect!!!');
  dum  = find(sum(W_data'));
  dum2 = find(sum(W_data));
  if isempty(A_first), A_firstred = []; else, A_firstred=A_first(dum,dum2); end
  if isempty(B_first), B_firstred = []; else, B_firstred=B_first(dum2,:); end
  [Ared,Bred,X_predred,info,A_list,B_list] = ...
      nca( X(dum,:),W_data(dum,dum2),W_signs(dum,dum2),'graphics_flag',graphics_flag,...
          'verbose',verbose,'force_B1',force_B1,'initial_choice',initial_choice,...
          'n_itmax',n_itmax,'lambda',lambda,'lambda_A',lambda_A,'lambda_B',lambda_B,...
          'A_first',A_firstred,'B_first',B_firstred,'repeat',repeat,'repeat_noise_std',repeat_noise_std,'epsilon',epsilon,...
          'no_nans',no_nans,'nan_blocks',nan_blocks,'TF_names',TF_names,...
          'general',general,'lambda_dep',lambda_dep,...
          'no_nan_predictions',no_nan_predictions,...
          'n_processes',n_processes,'dirname',dirname);
  
  [i,k] = size(W_data);
  l     = size(X,2);
  A = zeros(i,k);  A(dum,dum2) = Ared;
  B = randn(k,l);  B(dum2,:) = Bred;
  X_pred = zeros(i,l);  X_pred(dum,:) = X_predred;

  warning('Warning: irrelevant genes or TF. Output variables info, A_list, B_list are probably incorrect.');

else

  A_list = [];
  B_list = [];
  
  [i,k] = size(W_data);
  l     = size(X,2);

  random_seed = 1;
  [A_list,B_list] = nca_loop(X,W_data,W_signs,graphics_flag,verbose,force_B1,initial_choice,...
                             n_itmax,lambda,lambda_A,lambda_B,A_first,B_first,repeat,epsilon,...
                             no_nans,nan_blocks,TF_names,general,lambda_dep,regularise_flag,...
                             no_nan_predictions,n_processes,random_seed,repeat_noise_std);

end

% --- choose the n_accept best samples


% --- statistics over multiple runs

info = [];
info = nca_make_info(A_list,B_list,W_data,W_signs,X,lambda_A,force_B1,general,lambda_dep,regularise_flag);
B = info.B.best;
A = info.A.best;
B(B==0) = 10^-10;

A = optimise_A_given_B(X,W_data,B,W_signs,lambda_A,A,general,lambda_dep,regularise_flag);

X_pred = A*B;

if flag_X_positive,
  X_pred(find(X_pred<0))=0;
end

if no_nan_predictions, X_pred(isnan(X)) = nan; end

[dum,order] = sort(info.errors);
n_accept    = min(n_accept,length(order));
info.errors = info.errors(order(n_accept));
A_list      = A_list(order(1:n_accept));  
B_list      = B_list(order(1:n_accept));

end
