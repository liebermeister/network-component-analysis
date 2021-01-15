% [W,D] =  reorder_genes_W_D(gene_list,TF_list,W,D)

function [W,D] =  reorder_genes_W_D(gene_list,TF_list,W,D)

n_gene = length(W.gene_names);
n_TF     = length(W.TF_names);

if n_gene == n_TF, error('I am confused'); end

gene_order     = label_names(gene_list,W.gene_names,'single');
gene_inv_order = label_names(W.gene_names,gene_list,'single');
TF_order         = label_names(TF_list,W.TF_names,'single');
TF_inv_order     = label_names(W.TF_names,TF_list,'single');

gene_order = gene_order(find(gene_order));
TF_order = TF_order(find(TF_order));

% gene order

ff = fields(D);
for it = 1:length(ff),
  d = getfield(D,ff{it});
  if size(d) == n_gene,
    D=setfield(D,ff{it},d(gene_order));
  elseif size(d,1) == n_gene,
    D=setfield(D,ff{it},d(gene_order,:));
  end
end

ff = fields(W);
for it = 1:length(ff),
  d = getfield(W,ff{it});
  if size(d) == n_gene,
    W=setfield(W,ff{it},d(gene_order));
  elseif size(d,1) == n_gene,
    W=setfield(W,ff{it},d(gene_order,:));
  end
end

for it = 1:length(W.gene_sets),
  dum =  gene_inv_order(W.gene_sets{it});
  W.gene_sets{it} = sort(dum(find(dum)));
end

% TF order

ff = fields(W);
for it = 1:length(ff),
  d = getfield(W,ff{it});
  if length(d) == n_TF,
    W=setfield(W,ff{it},d(TF_order));
  elseif size(d,2) == n_TF,
    W=setfield(W,ff{it},d(:,TF_order));
  end
end

for it = 1:length(W.TF_sets),
  dum = TF_inv_order(W.TF_sets{it});
  W.TF_sets{it} = sort(dum(find(dum)));
end

dummy = find(cellfun('length',W.TF_sets) .* cellfun('length',W.TF_sets));
W.gene_sets = W.gene_sets(dummy);
W.TF_sets = W.TF_sets(dummy);