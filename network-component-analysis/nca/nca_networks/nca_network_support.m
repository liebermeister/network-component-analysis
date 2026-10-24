function M_support = nca_network_support(W,W_support)

% determine those connections in network W that are in addition supported by network W_support
ll_gene = label_names(W.gene_names,W_support.gene_names);
ll_TF   = label_names(W.TF_names,W_support.TF_names);
M_support = zeros(size(W.data));
M_support(find(ll_gene),find(ll_TF)) =  W.data(find(ll_gene),find(ll_TF)) .* W_support.data(ll_gene(find(ll_gene)),ll_TF(find(ll_TF)));
