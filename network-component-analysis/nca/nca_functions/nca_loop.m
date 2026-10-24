% [A_list,B_list] = nca_loop(...
%    X,W_data,W_signs,graphics_flag,verbose,force_B1,initial_choice,...
%    n_itmax,lambda,lambda_A,lambda_B,A_first,B_first,repeat,epsilon,...
%    no_nans,nan_blocks,TF_names,general,lambda_dep,regularise_flag,...
%    no_nan_predictions,n_processes,random_seed,repeat_noise_std)
%
% explanation of arguments: see nca.m

function [A_list,B_list] = nca_loop(...
    X,W_data,W_signs,graphics_flag,verbose,force_B1,initial_choice,...
    n_itmax,lambda,lambda_A,lambda_B,A_first,B_first,repeat,epsilon,...
    no_nans,nan_blocks,TF_names,general,lambda_dep,regularise_flag,...
    no_nan_predictions,n_processes,random_seed,repeat_noise_std)

rand( 'state',random_seed);
randn('state',random_seed);

if verbose, display([' Random seed: ' num2str(random_seed)]); end

for big_iteration = 1:repeat,
  
  if verbose, fprintf('Run %d/%d\n',big_iteration,repeat);  end
  
  iteration_finished = 0;
  
  while ~iteration_finished,

    try
  
      Xn = X + repeat_noise_std * randn(size(X));
  
      %% ----------------
      %% first guess 
      
      n_TF = size(W_data,2);
      
      clear A B goal_list resi_list;
      
      switch initial_choice,
        
        case 'given',
          
          A = A_first;
          B = B_first;
          
        case 'random',
          
          B = randn(size(W_data,2),size(X,2));
          A = W_data.*randn(size(W_data));
          A(find(A.*W_signs)<0) = - A(find(A.*W_signs)<0);
          
        case 'svd',
          
          for it=1:k,
            XX= Xn; XX(isnan(XX))=0;
            [u,s,v]= svd(XX(find(W_data(:,it)),:)); clear XX;
            A(find(W_data(:,it)),it) = u(:,1)*s(1,1);
            B(it,:)=v(:,1)';
            if sum(sign(B(it,:)))<0, B(it,:)=-B(it,:); A(:,it)=-A(:,it); end
          end
          
          d = sqrt(mean(B'.^2)).*sign(mean(B'))+10^-10;
          A = A * diag(d);
          B = diag(1./d) * B;
      
          B = B .* (1+0.5*randn(size(B))); % noise to resolve degeneracies
          if ~isempty(force_B1), B(1:size(force_B1,1),:)= force_B1; end
          B = diag(1./(10^-10+sqrt(mean(B'.^2)))) * B;
          
          if ~isempty(W_signs),  
            sgn = sign(sum(A.*W_signs)); 
            sgn(sgn==0)=1; A = A*diag(sgn); B = diag(sgn)*B; 
          end 
          
      end

      %% ----------------
      %% impose sign constraints
      
      ind_violate = find(A.*W_signs<0);
      if ind_violate,  A(ind_violate) = - 10^-3 * A(ind_violate);  end

      if graphics_flag, 
        figure(1); clf;
        subplot(1,2,1); im(B,[-2 2],TF_names);  title('B');
        drawnow;
        goal_list = [];
        resi_list = [];
        figure(2); clf;
      end
      
      %% ----------------
      %% optimisation
      
      if verbose, 
        fprintf('Updating matrices until rel. change < %f or %d iterations are reached\n',epsilon,n_itmax);
      end
      
      clear A_new B_new
      
      it   = 0;
      stop = (n_itmax==0);

      [err,err_residual, err_prior] = nca_goal(Xn,A,B,general,regularise_flag,lambda_dep,lambda_A,lambda_B);
      A_new = optimise_A_given_B(Xn,W_data,B,W_signs,lambda_A,A,general,lambda_dep,regularise_flag,lambda_B);
      
      while ~stop,
        
        it = it + 1;
        error_flag = 0;
        
        B_new = optimise_B_given_A(Xn,W_data,A_new,lambda_A,lambda_B,force_B1,no_nans,nan_blocks,...
                                   B,general,lambda_dep,regularise_flag);
        B_new = lambda * B_new + (1-lambda) * B;
        A_new = optimise_A_given_B(Xn,W_data,B_new,W_signs,lambda_A,A_new,general,lambda_dep,regularise_flag,lambda_B);
        A_new = lambda * A_new + (1-lambda) * A;
        [err_new, err_residual, err_prior] = nca_goal(Xn,A_new,B_new,general,regularise_flag,lambda_dep,lambda_A,lambda_B);
        err_change   = (err_new-err)/err;
        if err_new > err, error_flag = 1; end
        
        [A_new,B_new] = nca_rescale(A_new,B_new,force_B1);
        
        changes      = sqrt( mean(mean(  [B_new-B].^2 )) / mean(mean( B.^2 )) );
        
        if graphics_flag, 
          
          figure(1)
          subplot(1,2,1); im(B_new,[-2 2],TF_names);  title('B');
          subplot(1,2,2); im(B_new-B,[-0.2 0.2]);     title('Change of B');
          drawnow; 
          
          figure(2);
          goal_list = [goal_list err_new]; 
          resi_list = [resi_list err_residual]; 
          plot(goal_list); hold on; plot(resi_list,'r'); hold off; set(gca,'YScale','log');
          legend('Goal function','Residual');
          drawnow;
          
          figure(3); subplot(1,2,1); im(Xn); colorbar; subplot(1,2,2); im(A_new*B_new); colorbar;
          drawnow;
          
        end
        
        A    = A_new;
        B    = B_new;
        err  = err_new;
        
        if (abs(changes) < epsilon),  if abs(changes) > 10^-15, stop = 1;     end;   end
        
        %%    if ~stop, B = B + 10^-3/(it+1)^4 * randn(size(B)); end
        
        if verbose, 
          fprintf(' %d  Goal:%s  RelGoalChange:%s  Residual:%s  Prior term:%s RelChangeB:%s\n',...
                  it,num2str(err,3),num2str(err_change,3),num2str(err_residual,3),num2str(err_prior,3),num2str(changes,3));
          if error_flag,   fprintf('Error in NCA: goal function increases!\n');   end;
          if it == n_itmax, stop=1; 
            fprintf('No convergence reached. Increase maximal iteration number\n'); 
          end; 
        end
        
      end
      
      A     = optimise_A_given_B(Xn,W_data,B,W_signs,lambda_A,A,general,lambda_dep,regularise_flag,lambda_B);
      [A,B] = nca_rescale(A,B,force_B1);
      
      A_list{big_iteration} = A;
      B_list{big_iteration} = B;
      
      iteration_finished = 1; 
      
    catch
      display('nca_loop.m: An error occurred during optimisation; trying another starting point.');
    end
    
  end
  
end

if verbose,   fprintf(' Computation finished\n'); end
