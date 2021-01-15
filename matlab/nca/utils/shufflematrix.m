% function [smatrix] = shufflematrix(matrix,i)
% independently shuffle elements in columns (i=1) or in rows (i=2) 
% of a matrix

function [smatrix] = shufflematrix(matrix,i)

[m,n]=size(matrix);
smatrix = zeros(size(matrix));

if i==1,
  for k=1:n
   smatrix(:,k) = matrix(randperm(m),k);
  end
elseif i==2,
    for k=1:m
  smatrix(k,:) = matrix(k,randperm(n));
  end
end
