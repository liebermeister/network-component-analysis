function A_new = optimise_A_given_B(X,Wdata,B,signs,lambda_A,A,general,lambda_dep,regularise_flag,lambda_B)

% A_new = optimise_A_given_B(X,Wdata,B,signs,lambda_A,A,general,lambda_dep,regularise_flag)
%
%if ~exist('signs','var'),      signs       = []; end 
%if ~exist('general','var'),    general     = []; end 
%if ~exist('lambda_dep','var'), lambda_dep  = 1;  end 

ind_violate = find(A.*signs<0);
if ind_violate,  
  display('    Warning: Initial matrix A violates sign constraints'); 
  display('             Replacing all unfeasible entries by smaller feasible values'); 
  A(ind_violate) = - 10^-3 * A(ind_violate);
end

[i,k] = size(Wdata);
l     = size(B,2);

if regularise_flag,
  if length(lambda_dep) == 1,
    dumB = lambda_dep * B; 
    dumB(general,:) = 0;
  else,  
    dumB = diag(lambda_dep) * B;
  end
  X               = [X zeros(size(X))];
  B               = [B dumB];  
end

if lambda_A~=0,
  X = [X zeros(i,1)];
  B = [B lambda_A*ones(k,1)];
end

[i,k] = size(Wdata);

no_nans = sum(sum(isnan(X))) == 0;

if isempty(signs),  
  
  %% no sign restricion

  if no_nans,
    for it2=1:i,
      ind = find(Wdata(it2,:));
      A_new(it2,ind) = X(it2,:) * B(ind,:)'* pinv( B(ind,:) * B(ind,:)' + lambda_A );
    end
  else
    for it2=1:i,
      ind = find(Wdata(it2,:));
      valid = isfinite(X(it2,:));
      A_new(it2,ind) = X(it2,valid) * B(ind,valid)'* pinv( B(ind,valid) * B(ind,valid)' + lambda_A );
    end
  end
  
else,
  
  %% with sign restricion

  A_new = zeros(i,k);
  
  for it2=1:i,
    ind = find(Wdata(it2,:)); 
    s   = signs(it2,ind)';
    a_guess = A(it2,ind)';
    if no_nans,
      valid = 1:size(X,2);
      x = X(it2,:)';
      b = B(ind,:)';
    else,
      valid = isfinite(X(it2,:));
      x = X(it2,valid)';
      b = B(ind,valid)';
    end
    if find(s),
      A_new(it2,ind) = my_lsqlin(b,x,s,a_guess)';
    else, 
      A_new(it2,ind) = (pinv(b)*x)';
    end

    err1 = sum((    A(it2,ind) * b' - x').^2);
    err2 = sum((A_new(it2,ind) * b' - x').^2);
    if err2 > err1+10^-10, 
      fprintf('    Warning: A becomes worse\n'); 
      A_new(it2,ind) = A(it2,ind); 
    end
    
  end
  if find(sign(A_new).*signs<0),
    error;
    fprintf('Error in optimisation of A: sign constraint violated\n'); 
    A_new((A_new.*signs)<0)=0;
  end
  
end
