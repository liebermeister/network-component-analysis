function [Wsub,operon_indices,TF_indices,Wpicture] = pick_subW_operons(W,operon_list,sorting)

% [Wsub,operon_indices,TF_indices,Wpicture] = pick_subW_operons(W,operon_list,sorting)

if ~exist('sorting','var'), sorting = 1; end

operon_indices    = find_operon(operon_list,W);
TF_indices        = find(sum(W.data(operon_indices,:),1));

Wsub.data         = W.data(operon_indices,TF_indices);
Wsub.signs        = W.signs(operon_indices,TF_indices);
Wsub.gene_names   = W.gene_names(operon_indices);
Wsub.operon_names = W.operon_names(operon_indices);
Wsub.TF_names     = W.TF_names(TF_indices);

Wpicture          = Wsub.data;

if sorting, [Wsub,Wpicture] = sort_W_matrix(Wsub); end

%im(Wsub.data,[],Wsub.operon_names,Wsub.TF_names);
