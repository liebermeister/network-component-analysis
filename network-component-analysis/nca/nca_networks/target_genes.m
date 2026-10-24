function t = target_genes(TF,W)

% t = target_genes(TF,W)
% find the target genes of a TF
% examples: 
%  target_genes('argR',W)
%  target_genes({'argR','crp'},W)

if iscell(TF),
    t = {};
    for it =1:length(TF),
 t = [t; W.gene_names( find(W.data(:,find_TF(TF{it},W))))];
end
    t = unique(t);
else
t = W.gene_names( find(W.data(:,find_TF(TF,W))));
end
