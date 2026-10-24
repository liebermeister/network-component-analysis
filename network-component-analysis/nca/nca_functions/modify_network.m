function [A_opt,B_opt,W_opt,X_opt,X_error_opt,changes] = modify_network(A,B,W,X,X_error,method,parameters,general,lambda_dep,regularise_flag);

X_error(X_error==0)=1; %%%%%%%%%%%%%%% FIX IT !!!

A_opt = A;
B_opt = B;
W_opt = W;
X_opt = X; 
X_error_opt = X_error; 
changes = struct('sign_changes',[],'new',[],'removed',[]);  

for it2 = 1:parameters.iterations,
 for it = 1:length(method)
   switch method{it},
  
case 'refit_A',
  A_opt = optimise_A_given_B(X_opt,W_opt.data,B_opt,W_opt.signs,0);

case 'shift_time',
 [X_opt,changes.shifts,A_opt,X_error_opt] = modify_shifts(W_opt,X_opt,A_opt,B_opt,10,[],X_error_opt,parameters.experiments);
% if parameters.verbose, fprintf('Shifted:\n');   end
 
case 'relax_sign',
  [W_opt,A_opt,changes.sign_changes] = modify_signs(X_opt,W_opt,A_opt,B_opt,W_opt.signs,X_error_opt,parameters.p_sign_error,parameters.fixed_TF,0,general,lambda_dep,regularise_flag);

case  'change_edges'
 [W_opt,A_opt,X_pred,B,changes.new,changes.removed] = optimise_W(X_opt,parameters.W_match,A_opt,B_opt,parameters.P,parameters.W_match.signs,parameters.sigma_0,parameters.p_sign_error,parameters.n_threshold,1,0,0,general,lambda_dep,regularise_flag);

  end
end
end

[A_opt,B_opt,X_opt] = nca(X,W_opt.data,W_opt.signs,'graphics_flag',1,'nan_blocks',parameters.experiments,'epsilon',0.0001,'TF_names',W.TF_names,'initial_choice','given','A_first',A_opt,'B_first',B_opt);

for it = 1:length(changes.new),
  changes.new_sign(it,1) = sign(A_opt(changes.new(it,2),changes.new(it,1) ));
end
