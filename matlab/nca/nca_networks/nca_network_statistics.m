function nca_network_statistics(W, B, filename)

ind_TF_fitted = find(sum(B>10^-5,2))';

eval(default('filename','[]'));
display(sprintf('Network statistics: %d genes, %d transcription factors, %d interactions, %d activations, %d repressions, %d interactions of unknown sign\n', length(W.gene_names), length(W.TF_names(ind_TF_fitted)), sum(sum(W.data(:,ind_TF_fitted)~=0)), sum(sum(W.signs(:,ind_TF_fitted)==1)), sum(sum(W.signs(:,ind_TF_fitted)==-1)), sum(sum((W.data(:,ind_TF_fitted)~=0).*(abs(W.signs(:,ind_TF_fitted))~=1)))))

if length(filename),

  T = {{'Element','Count number'};...
       {'Genes',                       num2str(length(W.gene_names))};...
       {'Transcription factor',        num2str(length(W.TF_names))};...
       {'Transcription factor varying',num2str(length(W.TF_names(ind_TF_fitted)))};...
       {'Interaction',                 num2str(sum(sum(W.data(:,ind_TF_fitted)~=0)))};...
       {'Activation',                  num2str(sum(sum(W.signs(:,ind_TF_fitted)==1)))};...
       {'Repression',                  num2str(sum(sum(W.signs(:,ind_TF_fitted)==-1)))};...
       {'Interactions of unknown sign',num2str(sum(sum((W.data(:,ind_TF_fitted)~=0).*(abs(W.signs(:,ind_TF_fitted))~=1))))}};
  
 mytable(T,0,filename);
  display(sprintf('Writing file %s',filename));  
end