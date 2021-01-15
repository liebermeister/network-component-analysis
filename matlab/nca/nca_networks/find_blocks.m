function [TF_sets,operon_sets] = find_blocks(W)

dist = W.data * W.data' ;
g.matrix = double(dist>0);
[operon_sets,node_order] = find_separate_subgraphs(g);

for it =1:length(operon_sets)
  TF_sets{it} = find(sum(W.data(operon_sets{it},:),1));
end