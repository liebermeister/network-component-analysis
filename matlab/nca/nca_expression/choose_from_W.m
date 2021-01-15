% W = choose_from_W(W,l_gene,lTF)
function W = choose_from_W(W,l_gene,lTF)


W.TF_names     = W.TF_names(lTF);
W.data         = W.data(:,lTF);
W.signs        = W.signs(:,lTF);
W.gene_names   = W.gene_names(l_gene);
W.operon_names = W.operon_names(l_gene);
W.operon_flags = W.operon_flags(l_gene);
W.operon_abbr  = W.operon_abbr(l_gene);
W.data         = W.data(l_gene,:);
W.signs        = W.signs(l_gene,:);
if isfield(W,'gene_sets'),
  W              = rmfield(W,'gene_sets');
  W              = rmfield(W,'TF_sets');
end

