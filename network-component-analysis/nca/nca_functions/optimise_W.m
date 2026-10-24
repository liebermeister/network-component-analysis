function [W_new,A_new,X_pred,B,new,removed] = optimise_W(X,W_start,A,B,P,known_signs,sigma_0,p_sign_error,n_threshold_possible,n_iterations,graphics_flag,verbose,general,lambda_dep,regularise_flag);

% [W_new,A_new,X_pred,B,new,removed] = optimise_W(X,W_start,A,B,P,known_signs,sigma_0,p_sign_error,n_threshold_possible,n_iterations,graphics_flag,verbose);
%
% Iterate over the adaptation of W_new and A,B

lambda_regular_A = 0;
if ~exist('graphics_flag','var'), graphics_flag = 0; end 

W_new = W_start;

if graphics_flag, 
 figure(1); 
 subplot(1,3,1); im(A*B); 
 subplot(1,3,2); im(A); 
 subplot(1,3,3); im(B); drawnow
end

for it =1:n_iterations,

  fprintf('Updating the connectivity matrix: iteration %d\n',it);
  
  W_this = modify_connections(X,W_new,A,B,P,known_signs,sigma_0,p_sign_error,n_threshold_possible,...
                              verbose,general,lambda_dep,regularise_flag);
  
  if graphics_flag,
    figure(2)
    subplot(2,2,1); im(W_start.data); title('START');
    subplot(2,2,2); im(W_new.data);       title('OLD');
    subplot(2,2,3); im(W_this.data);   title('UPDATED');
    drawnow
  end

  W_new = W_this;
  if it < n_iterations, [A,B] = nca(X,W.data,W.signs); end
  
  if graphics_flag, 
    figure(3); subplot(1,3,1); im(A*B); subplot(1,3,2); im(A); subplot(1,3,3); im(B); drawnow; 
  end
  
end

A_new = optimise_A_given_B(X,W_new.data,B,W_new.signs,lambda_regular_A,A,general,lambda_dep,regularise_flag);
if nargout > 3, X_pred = A_new*B; end

shifts = [];

dummy = (W_new.data - W_start.data );
 [i1,i2] = ind2sub(size(W_start.data),find(dummy>0));
 new = [i2 i1];
 [i1,i2] = ind2sub(size(W_start.data),find(dummy<0));
 removed = [i2,i1];
 
 if verbose,
   fprintf('%d new connections, %d connections removed\n',sum(dummy(:)>0),sum(dummy(:)<0));
   'new'
   [W_start.TF_names(new(:,1)) W_start.gene_names(new(:,2)) ]
   'removed'
   [W_start.TF_names(removed(:,1)) W_start.gene_names(removed(:,2)) ]
 end
 