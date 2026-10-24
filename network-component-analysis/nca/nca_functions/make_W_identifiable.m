function [Wblock,genes_kept,TF_removed,genes_removed] =  make_W_identifiable(W,gene_indices_choice,TF_indices_choice, general)

% [Wblock,genes_kept,TF_removed,genes_removed] =  make_W_identifiable(W,gene_indices_choice,TF_indices_choice, general)
%
% extract identifiable TF subnetwork from given network
% only argument W is required

if ~exist('general','var'), general = 0; end
if ~exist('gene_indices_choice','var'), gene_indices_choice = 1:size(W.data,1); end
if ~exist('TF_indices_choice','var'), TF_indices_choice = 1:size(W.data,2); end

Wblock.data         = W.data(gene_indices_choice,TF_indices_choice);
Wblock.signs        = W.signs(gene_indices_choice,TF_indices_choice);
Wblock.TF_names     = W.TF_names(TF_indices_choice);
Wblock.gene_names   = W.gene_names(gene_indices_choice);
Wblock.operon_names = W.operon_names(gene_indices_choice);
Wblock.operon_abbr  = W.operon_abbr(gene_indices_choice);
Wblock.operon_flags = W.operon_flags(gene_indices_choice);

% --- check rank conditions

er = 1;

genes_kept  = 1:size(Wblock.data,1);
genes_removed = [];
TF_removed    = [];

[er,independent_TF,dependent_TF,T,violating_columns,dependent_columns] = nca_identifiability(Wblock.data.*randn(size(Wblock.data)));

while er,

dependent_operons   = find(sum(Wblock.data(:,dependent_TF),2));
independent_operons = setdiff(1:size(Wblock.data,1),dependent_operons);

TF_removed      = [TF_removed; Wblock.TF_names(dependent_TF)];
genes_removed   = [genes_removed; Wblock.gene_names(dependent_operons)];
genes_kept      = genes_kept(independent_operons); 

Wblock.TF_names     = Wblock.TF_names(independent_TF);
Wblock.gene_names   = Wblock.gene_names(independent_operons);
Wblock.operon_names = Wblock.operon_names(independent_operons);
Wblock.operon_abbr  = Wblock.operon_abbr(independent_operons);
Wblock.operon_flags = Wblock.operon_flags(independent_operons);
Wblock.data         = Wblock.data(independent_operons,independent_TF);
Wblock.signs        = Wblock.signs(independent_operons,independent_TF);

[er,independent_TF,dependent_TF,T,violating_columns,dependent_columns] = nca_identifiability(Wblock.data.*randn(size(Wblock.data)));

independent_columns = setdiff(1:size(Wblock.data,2), dependent_columns);
keep_genes        = setdiff(1:size(Wblock.data,1), find(sum(Wblock.data(:,dependent_columns)')));
genes_kept = genes_kept(keep_genes);
fprintf(' Matrices not yet identifiable. Searching identifiable subnetwork.\n');
TF_removed = [TF_removed; Wblock.TF_names(dependent_columns)];
genes_removed   = [genes_removed; Wblock.gene_names(find(sum(Wblock.data(:,independent_columns)')) )];
Wblock.TF_names = Wblock.TF_names(independent_columns);
Wblock.data = Wblock.data(keep_genes,independent_columns);
Wblock.signs = Wblock.signs(keep_genes,independent_columns);
Wblock.operon_names = Wblock.operon_names(keep_genes);
Wblock.operon_abbr  = Wblock.operon_abbr(keep_genes);
Wblock.operon_flags = Wblock.operon_flags(keep_genes);
Wblock.gene_names   = Wblock.gene_names(keep_genes);

[er,independent_TF,dependent_TF,T,violating_columns,dependent_columns] = nca_identifiability(Wblock.data.*randn(size(Wblock.data)));

end

% --- add general input

if general,
if max(sum(Wblock.data))<size(Wblock.data,1),
Wblock.TF_names =[{'General input'}; Wblock.TF_names];
Wblock.data = [ones(size(Wblock.data,1),1) Wblock.data];
Wblock.signs = [ones(size(Wblock.signs,1),1) Wblock.signs];
end
end
