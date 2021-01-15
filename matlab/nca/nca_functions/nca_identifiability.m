function [error_flag,independent_TF,dependent_TF,T,violating_columns,dependent_columns] = nca_identifiability(A_true,B_true)

% [error_flag,independent_TF,dependent_TF,T,violating_columns,dependent_columns] = nca_identifiability(A_true,B_true)
% error_flag      1 if matrices are unidentifiable
%
% first condition: 
%   independent_TF  indices of TF not violating the first condition
%   dependend_TF
%   T

A_true = full(A_true);

k = size(A_true,2);

error_flag = 0;

violating_columns = [];

% condition 1

dependent_TF = [];
independent_TF = 1:size(A_true,2);
T=eye(size(A_true,2));
if k ~= min(rank(A_true),k), 
  fprintf(' Condition 1 violated. Computing set of independent columns\n'); error_flag = 1; 
  [echelon,independent_TF] = rref(A_true);
  T                    = echelon(1:length(independent_TF),:);
dependent_TF = setdiff(1:size(A_true,2),independent_TF);
end

% condition 2

  violating_columns = [];
  dependent_columns = [];
 
if error_flag==0,
% condition 2
 for it =1:k,
   dum=A_true(find(A_true(:,it)==0),setdiff(1:k,it));
     if k-1 ~= min(rank(dum),k-1), fprintf(' Condition 2 violated by column %d\n',it); error_flag = 1;
          violating_columns=[  violating_columns; it];
          [echelon,independent_col] = rref(dum);
          dependent_col = setdiff(1:size(dum,2),independent_col);
	   column_choice = setdiff(1:k,it);
           dependent_columns = [dependent_columns column_choice(dependent_col)];
     end
  end
end

dependent_columns = unique(dependent_columns);

% condition 3
if exist('B','var'),
if k ~= min(rank(B_true),k), fprintf(' Condition 3 violated\n'); error_flag = 1; end
end
  if error_flag==0, fprintf(' The matrices are identifiable\n');
end
