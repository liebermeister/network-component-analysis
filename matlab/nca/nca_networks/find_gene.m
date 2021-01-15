function index = find_gene(gene,D)

if isstr(gene),
index = label_names({gene},D.gene_names,'multiple');
index = index{1};
else,
index = label_names(gene,D.gene_names);
end

if isempty(index), fprintf('gene %s not found\n',gene); end