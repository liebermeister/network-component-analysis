function [Wsub,gene_indices,TF_indices,Wpicture] = pick_subW_genes(W,gene_list,sorting)

% [Wsub,gene_indices,TF_indices,Wpicture] = pick_subW_genes(W,gene_list,sorting)

if ~exist('sorting','var'), sorting = 1; end

gene_indices      = find_gene(gene_list,W);
gene_indices      = gene_indices(gene_indices~=0);

TF_indices        = find(sum(W.data(gene_indices,:),1));

Wsub.data         = W.data(gene_indices,TF_indices);
Wsub.signs        = W.signs(gene_indices,TF_indices);
Wsub.gene_names   = W.gene_names(gene_indices);
Wsub.operon_names = W.operon_names(gene_indices);
Wsub.TF_names     = W.TF_names(TF_indices);

Wpicture          = Wsub.data;

if sorting, [Wsub,Wpicture] = sort_W_matrix(Wsub); end

%im(Wsub.data,[],Wsub.gene_names,Wsub.TF_names);
