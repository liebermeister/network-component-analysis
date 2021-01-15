function [p,q] = sort_block(M,only_rows,nit)

eval(default('only_rows','0','nit','[]'));

% reduce the problem by replacing groups of identical rows by a single row each

MM = unique(M,'rows'); 

if size(MM,1)<size(M,1),
  [p,q] = sort_block(MM,only_rows,nit);
  pp = [];
%   z  = 0;
  for it = 1:length(p)
    ind = find(sum(abs(repmat(MM(p(it),:),size(M,1),1)-M),2)==0);
    pp = [pp; ind];
%     z = z+length(ind);
  end
  p = pp;
else,

  if isempty(nit),
    nit = ceil(100+5*prod(size(M)));
  end
  
  display(sprintf('Solving traveling salesman problem, %d iterations',nit)); 
  
  dist_OP = M*M';
  seed    = 0;
  
  if size(M,1)>1,
    p = traveling_salesman(M'*diag(1./sqrt(diag(dist_OP))),nit,'euclidean',seed);
    if (size(M,2)>2) * (only_rows==0),
      M = M(p,:);
      dist_TF = M'*M;
      q = traveling_salesman(M*diag(1./sqrt(diag(dist_TF))),nit,'euclidean');
    else,  
      q = 1:size(M,2); 
    end
  else,
    p = 1;
    q = 1:size(M,2); 
  end
  
  
  %% sort columns by their center of gravity (along the column)
  score = mean([M~=0] .* repmat([1:size(M,1)]',1,size(M,2))) ./ sum(M~=0);
  [dum, q] = sort(score);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

M = M(p,q);

M_old = M;
stop = 0;

while ~stop,
  M_old = M;
  
  A   = repmat(1:size(M,2),size(M,1),1);
  dum = mean(diag(1./sum(M,2))*A.*M,2);
  [dum2,order] = sort(dum);
  p = p(order);
  M = M(order,:);
  
  if ~only_rows,
    A   = repmat((1:size(M,1))',1,size(M,2));
    dum = mean(A.*M*diag(1./sum(M,1)));
    [dum2,order]=sort(dum); 
    q = q(order);
    M = M(:,order);
  end
  
  if sum(sum(abs(M-M_old)))==0, stop = 1; end
end

%end