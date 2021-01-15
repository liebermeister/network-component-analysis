function [subgraphs,node_order] = find_separate_subgraphs(g)

g_size=size(g.matrix,1);
M = expm(full(g.matrix))>0;

nl=1:g_size;

subgraphs={};
z=1;
while ~isempty(nl),
       pattern=M(:,nl(1));
       subgraphs{z}= find(sum(abs(M-repmat(pattern,1,g_size)))==0);
       nl=setdiff(nl, subgraphs{z});
       z=z+1;
end

[dum,order]=   sort( -cellfun('length',subgraphs));

subgraphs=subgraphs(order);
     
node_order = [];
for it =1:length(subgraphs), node_order = [node_order subgraphs{it}]; end 
