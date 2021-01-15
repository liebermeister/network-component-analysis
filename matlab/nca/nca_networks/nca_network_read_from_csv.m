function W = nca_network_read_from_csv(filename)

M = load_any_table(filename);

W.gene_names = unique(M(2:end,1));
W.TF_names   = unique(M(2:end,2));

ind_gene = label_names(M(2:end,1),W.gene_names);
ind_TF   = label_names(M(2:end,2),W.TF_names);

W.data  = sparse(zeros( length(W.gene_names), length(W.TF_names) ));
W.signs = W.data;

ind_all = sub2ind(size(W.data),ind_gene,ind_TF);
W.data(ind_all) = 1;
W.signs(ind_all(strcmp(M(2:end,5),'+1'))) = 1; 
W.signs(ind_all(strcmp(M(2:end,5),'-1'))) = -1; 

W.operon_names = W.gene_names; 
% Fake entries since actual operon names are not known