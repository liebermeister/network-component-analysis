function W_data_count = W_statistics(W,Wsparse_list,Wsparse_names)

% statistics over a network structure and some filtered versions of it

display(sprintf('Entire network:'));
display(sprintf(' Number of TF   : %d',length(W.TF_names)));
display(sprintf(' Number of genes: %d',length(W.gene_names)));
display(sprintf(' Number of edges: %d',sum(W.data(:))));

if exist('Wsparse_list','var'),
  W_data_count     = zeros(size(W.data));
  for it = 1:length(Wsparse_list),
    display(sprintf('Filtered network %s: %d edges (fraction: %0.3f)', Wsparse_names{it}, sum(Wsparse_list{it}.data(:)), sum(Wsparse_list{it}.data(:)) /sum(W.data(:)) ));
    W_data_count = W_data_count + Wsparse_list{it}.data;
  end
  W_data_intersect = W_data_count == length(Wsparse_list);
  W_data_consensus = W_data_count>=ceil([length(Wsparse_list)/2]);
  W_data_union     = W_data_count>0;
  display(sprintf('Intersection network: %d edges (fraction: %0.3f)',sum(W_data_intersect(:)) ,sum(W_data_intersect(:)) /sum(W.data(:)) ));
  display(sprintf('Union network       : %d edges (fraction: %0.3f)',sum(W_data_union(:))     ,sum(W_data_union(:)) /sum(W.data(:)) ));
  display(sprintf('Consensus network   : %d edges (fraction: %0.3f)',sum(W_data_consensus(:)),sum(W_data_consensus(:)) /sum(W.data(:)) ));
  display(sprintf(' ("Consensus" means: edges present in at least %d of the filtered networks)',ceil([length(Wsparse_list)/2])));
end