function [emat,evec] = make_emat(nan_blocks,indices,predict_matrix,Xtrue,Xpred)

for it =1:length(nan_blocks),
  for it2=1:length(indices),
    emat(it2,it) = nanmean( ((Xtrue(indices(it2),nan_blocks{it})-Xpred(indices(it2),nan_blocks{it})).^2)' );
  end 
  emat(predict_matrix(indices,nan_blocks{it}(1))==0,it)=nan;
end

for it2=1:length(indices),
  ind_i = indices(it2);
  ind_k = find(predict_matrix(indices(it2),:));
  evec(it2,1) =  nanmean( (Xtrue(ind_i,ind_k) - Xpred(ind_i,ind_k)).^2' );
end
