% [W,D] =  reorder_genes(operon_list,TF_list,W,D)

function [A,B,W,X_pred,X,Xplus,Xminus,X_error,info] =  reorder_genes(operon_list,TF_list,X,Xmean,X_error,Xplus,Xminus,X_pred,W,A,B,info,D)

gene_order = label_names(operon_list,W.operon_names,'single');
TF_order   = label_names(TF_list,W.TF_names,'single');

if exist('A','var'), A=A(gene_order,:); end

ff = fields(D);
for it = 1:length(ff),
  d = getfield(D,ff{it});
  if size(d) == length(gene_order),
    D=setfield(D,ff{it},d(gene_order));
  elseif size(d,1) == length(gene_order),
    D=setfield(D,ff{it},d(gene_order,:));
  end
end

ff = fields(W);
for it = 1:length(ff),
  d = getfield(W,ff{it});
  if size(d) == length(gene_order),
    W=setfield(W,ff{it},d(gene_order));
  elseif size(d,1) == length(gene_order),
    W=setfield(W,ff{it},d(gene_order,:));
  end
end

for it = 1:length(W.operon_sets),
  W.operon_sets{it} = gene_order(W.operon_sets{it});
end

X       =       X(gene_order,:);
Xplus   =   Xplus(gene_order,:);
Xminus  =  Xminus(gene_order,:);
X_error = X_error(gene_order,:);
X_pred  =  X_pred(gene_order,:);
Xmean   =   Xmean(gene_order);

if exist('info','var'),
  iinfo.B = info.B;
  iinfo.errors = info.errors;
  iinfo.A.mean = info.A.mean(gene_order,:);
  iinfo.A.std  = info.A.std(gene_order,:);
  iinfo.A.best = info.A.best(gene_order,:);
  for it=1:length(info.A.list), iinfo.A.list{it} = info.A.list{it}(gene_order,:); end
info = iinfo;
end

clear operon_indices_choice

% TF order

if exist('A','var'), A=A(:,TF_order); end
B=B(TF_order,:);

ff = fields(W);
for it = 1:length(ff),
  d = getfield(W,ff{it});
  if length(d) == length(TF_order),
    W=setfield(W,ff{it},d(TF_order));
  elseif size(d,2) == length(TF_order),
    W=setfield(W,ff{it},d(:,TF_order));
  end
end

for it = 1:length(W.TF_sets),
  W.TF_sets{it} = TF_order(W.TF_sets{it});
end

if exist('info','var'), 
  iinfo.errors = info.errors;
  iinfo.A.mean = info.A.mean(:,TF_order);
  iinfo.A.std  = info.A.std(:,TF_order);
  iinfo.A.best = info.A.best(:,TF_order);
  iinfo.B.mean = info.B.mean(TF_order,:);
  iinfo.B.std  = info.B.std(TF_order,:);
  iinfo.B.best = info.B.best(TF_order,:);
  for it=1:length(info.A.list), iinfo.A.list{it} = info.A.list{it}(:,TF_order); end; 
  for it=1:length(info.B.list), iinfo.B.list{it} = info.B.list{it}(TF_order,:); end
info = iinfo;
end
