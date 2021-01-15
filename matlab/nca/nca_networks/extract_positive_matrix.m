function [keep_operons,keep_TF] = extract_positive_matrix(W);

% [keep_operons,keep_TF] = extract_positive_matrix(W);
% extract a subnetwork with nonnegative connections
% (useful for non-negative matrix factorisation)

A = zeros(size(W.signs));
A(W.signs==1) = 1;
A(W.signs>1) = -1;

%im(A)

keep_operons = 1:size(A,1);
keep_TF = 1:size(A,1);

col_score=sum(A<0)./(1+sum(A>0));
row_score=sum(A'<0)./(1+sum(A'>0));

threshold = 1;
while max(col_score)>0 | max(row_score)>0,
 threshold = threshold * 0.95;
 k_TF = find(col_score<=nanquantile(col_score',threshold));
 r_TF = setdiff(1:size(A,2),k_TF);
 k_operons = find(sum(abs(A(:,r_TF))')==0);
 A=A(k_operons,k_TF);

 keep_operons = keep_operons(k_operons);
 keep_TF =keep_TF(k_TF);

 col_score=sum(A<0)./(1+sum(A>0));
 row_score=sum(A'<0)./(1+sum(A'>0));

 A=A(find(row_score<=nanquantile(row_score',threshold)),:);

keep_operons = keep_operons(find(row_score<=nanquantile(row_score',threshold)));

col_score=sum(A<0)./(1+sum(A>0));
row_score=sum(A'<0)./(1+sum(A'>0));

%im(A)
% pause;
end

keep_operons = keep_operons(find(sum(A')>0));
keep_TF = keep_TF(find(sum(A)>0));
