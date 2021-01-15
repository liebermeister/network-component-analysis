function [W, gene_order, TF_order, Wpicture] = sort_W_matrix(W,only_rows)


% [W,Wpicture,gene_order,TF_order] = sort_W_matrix(W)
% reorder genes and TF in order to get a nice "almost blockdiagonal" matrix

% disregard first column if it is full

eval(default('only_rows','0'));

if (sum(W.data(:,1)==0)==0) * ( size(W.data,2)>1),

  dumW          = W;
  dumW.data     = W.data(:,2:end);
  dumW.signs    = W.signs(:,2:end);
  dumW.TF_names = W.TF_names(2:end);

  [dumW, gene_order, TF_order, dumWpicture] = sort_W_matrix(dumW);
  Wpicture = [2*ones(length(W.operon_names),1) dumWpicture ];

  TF_order = [1, 1+column(TF_order)'];
  W.data   = [W.data(:,1) dumW.data];
  W.signs  = [W.data(:,1) dumW.signs];
  W.operon_names = dumW.operon_names;

  if isfield(W,'gene_names'),    W.gene_names   = dumW.gene_names;   end
  if isfield(W,'operon_abbr'),   W.operon_abbr  = dumW.operon_abbr;  end
  if isfield(W,'operon_flags'),  W.operon_flags = dumW.operon_flags; end

  W.TF_names  = [W.TF_names(1); dumW.TF_names];
  W.gene_sets = dumW.gene_sets;
  W.TF_sets   = [{1} dumW.TF_sets];

  for it =1:length(W.TF_sets), W.TF_sets{it} = W.TF_sets{it}+1; end
  
else,

%   if only_rows,
%     TF_order = 1:length(W.TF_names);
%     nit = ceil(100+5*prod(size(W.data)));
%     gene_order =  sort_block(W.data,1,nit);
%     operon_set_size        = length(W.gene_names);
%     TF_set_size            = length(W.TF_names);
%   else,
    
% --- define blocks in W

% -- indices in lists operon_set and TF_set

  dist = W.data * W.data' ;
  g.matrix = double(dist>0);
  [subgraphs,node_order] = find_separate_subgraphs(g);
  
  operon_set = subgraphs;
  TF_set={};
  for it = 1:length(subgraphs), 
    TF_set{it} = find(sum(W.data(subgraphs{it},:),1)>0);
  end 
  
  [dum,order] = sort(-cellfun('length',TF_set));
  operon_set  = operon_set(order);
  TF_set      = TF_set(order);
  
% --- order within each block
  
  for it = 1:length(operon_set),
    if length(TF_set{it}) > 1 & length(operon_set{it}) > 1  ,
      block          = W.data(operon_set{it},TF_set{it});
%      size(block)
%      if size(block,2)>1,
        [p,q]          = sort_block(block,only_rows);
%      else,
%        p = [1:size(block,1)]'; q = 1;
%      end
      operon_set{it} = operon_set{it}(p);
      TF_set{it}     = TF_set{it}(q); 
    end
end

% --- reorder W according to blocks

TF_order   = [];
gene_order = [];

for it =1:length(subgraphs), 
  TF_order   = [TF_order, column(TF_set{it})'];
  gene_order = [ gene_order column(operon_set{it})' ];
  TF_set_size(it)     = length(TF_set{it});
  operon_set_size(it) = length(operon_set{it});
end 

if only_rows, TF_order = 1:length(TF_order); end

W.data         = W.data(gene_order,TF_order);
W.signs        = W.signs(gene_order,TF_order);
W.operon_names = W.operon_names(gene_order);
if isfield(W,'gene_names'),   W.gene_names   = W.gene_names(gene_order); end
if isfield(W,'operon_abbr'),  W.operon_abbr  = W.operon_abbr(gene_order); end
if isfield(W,'operon_flags'), W.operon_flags = W.operon_flags(gene_order); end

W.TF_names = W.TF_names(TF_order);

% show block matrix

Wpicture = W.data;
s1 = cumsum([0 operon_set_size]);
s2 = cumsum([0 TF_set_size]);
for it = 1:length(s1)-1;
  operon_indices{it} = s1(it)+1:s1(it+1);
  TF_indices{it} = s2(it)+1:s2(it+1);
  Wpicture(operon_indices{it},TF_indices{it}) = Wpicture(operon_indices{it},TF_indices{it})+1;
end

W.gene_sets = operon_indices;
W.TF_sets   = TF_indices;

end

