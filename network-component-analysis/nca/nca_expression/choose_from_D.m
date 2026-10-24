% D = choose_from_D(D,indices)

function D = choose_from_D(D,indices)

n_gene = length(D.gene_names);

ff = fields(D);
for it = 1:length(ff),
  d = getfield(D,ff{it});
  if size(d) == n_gene,
    D=setfield(D,ff{it},d(indices));
  elseif size(d,1) == n_gene,
    D=setfield(D,ff{it},d(indices,:));
  end
end
