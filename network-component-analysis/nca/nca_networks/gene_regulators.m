function t = gene_regulators(gene,W)

% t = regulators(gene,W)
% find the regulators of a gene or a list of genes

if iscell(gene),
  t = {};
  for it = 1:length(gene),
    this_t = gene_regulators(gene{it},W);
    t = [t; this_t];
  end
  t = unique(t);
else
t = W.TF_names( find(W.data(find_gene(gene,W),:)));

end
