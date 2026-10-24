function [D_comb,W_comb,W_picture,D_unmatched] = combine_D_and_W(D,W)

% [D,W,W_picture,D_unmatched] = combine_D_and_W(D,W)
%
% Combine experimental data in D with network in W
% according to the GENE NAMES!!!

eval('!rm /tmp/map_geneD_geneW.txt');
eval('!rm /tmp/Wnames.txt');
eval('!rm /tmp/Dnames.txt');

all_genes = W.gene_names;
table(all_genes,1,'/tmp/Wnames.txt');

all_genes = D.gene_names;
table(all_genes,1,'/tmp/Dnames.txt');

eval('! ~/projekte/NCA/Experimental_data/Expression_data/Constructs/match_names_alon_ecocyc.pl > /tmp/map_geneD_geneW.txt');
map_geneD_geneW = load('/tmp/map_geneD_geneW.txt');

%if isfield(D,'allgene_names'), 
%  dum = zeros(length(D.gene_names),2);
%  dum(:,1)=(1:length(D.gene_names))';
%  for it=1:length(D.gene_names),
%    ll = label_names(D.operon2gene{it},D.allgene_names);
%    res = map_geneD_geneW(ll,2);
%    res = res(res>0);
%    if res,
%      dum(it,2) = res(1);
%    end 
%  end
%  map_geneD_geneW = dum;
%end

Dlines = find(map_geneD_geneW(:,2));

if ~isfield(D,'operon_names'), 
  [D.operon_names,D.operon_flags,D.operon_abbr] = e_coli_gene2operon(D.gene_names);
end 

%[D.operon_names(Dlines) D.gene_names(Dlines) W.gene_names(map_geneD_geneW(Dlines,2))];

unmatched_genes = find(map_geneD_geneW(:,2)==0);
D_unmatched = choose_from_D(D,unmatched_genes);

D.map_geneD_geneW = map_geneD_geneW(:,2);

% reduce W matrix to available expression data

lines = D.map_geneD_geneW(find(D.map_geneD_geneW));

Wred = W;
Wred.gene_names   = D.gene_names(Dlines);
Wred.operon_names = D.operon_names(Dlines);
Wred.operon_flags = D.operon_flags(Dlines);
Wred.operon_abbr  = D.operon_abbr(Dlines);
Wred.data         = Wred.data(lines,:);
Wred.signs        = Wred.signs(lines,:);

% keep only relevant TF

dum = find(sum(Wred.data));
Wred.data = Wred.data(:,dum);
Wred.signs = Wred.signs(:,dum);
Wred.TF_names = Wred.TF_names(dum);

%figure(1); 
%im(Wred.data,[],Wred.gene_names,Wred.TF_names)
%set(gca,'FontSize',8)

% -----------------------------------
% reordering W

[Wsorted,ind_operon,ind_TF,W_picture] = sort_W_matrix(Wred);

% choose only operons that are present in shalevs network
dum   = find(D.map_geneD_geneW); 
order = dum(ind_operon);

Dsorted   = choose_from_D(D,order);

if isfield(D,'original_indices'),
  Dsorted.map_geneD_geneW = D.map_geneD_geneW(order); 
end

%figure(2); im(W_picture,[],Dsorted.gene_names,Wsorted.TF_names);
%set(gca,'FontSize',8)

% ----------------------------------

%Wsorted.gene_names = Dsorted.gene_names;

W_comb = Wsorted;
D_comb = Dsorted;
