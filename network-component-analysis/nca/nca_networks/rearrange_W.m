function W_new = rearrange_W(W,gene_names,TF_names)

% W_new = rearrange_W(W,gene_names,TF_names)
%
% reorder network / pick subnetwork
% according to lists of genes and TF
%
% if genes or TF are not found, then
% empty rows and columns are created

l1 = label_names(gene_names,W.gene_names);
l2 = label_names(TF_names,W.TF_names);

W_new.gene_names               = cell(length(gene_names),1);
W_new.operon_names             = cell(length(gene_names),1);
W_new.TF_names                 = cell(length(TF_names),1);
W_new.gene_names               = gene_names;
W_new.operon_names(find(l1))   = W.operon_names(l1(find(l1)));
W_new.TF_names(find(l2))       = TF_names(l2(find(l2)));
W_new.data = zeros(length(gene_names),length(TF_names));
W_new.data(find(l1),find(l2))  = W.data(l1(find(l1)),l2(find(l2)));
W_new.signs = zeros(length(gene_names),length(TF_names));
W_new.signs(find(l1),find(l2)) = W.signs(l1(find(l1)),l2(find(l2)));

