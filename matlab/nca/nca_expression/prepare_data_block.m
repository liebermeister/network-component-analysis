function [D,W,TF_removed,gene_indices_choice] = prepare_data_block(filename,identifiable_flag,TF_list,gene_list)

% prepare_data_block(filename,identifiable_flag,TF_list,gene_list)
%
% choose a part of the whole data, display and save them

cd_data(filename,'analysis');
load(['matched_' filename]);


% ---------------------------------

%figure(1); im(W_picture,[],W.gene_names,W.TF_names); set(gca,'Fontsize',6);

% --------------------------------
% choose all data ...

gene_indices_choice = 1:length(W.gene_names);
TF_indices_choice     = 1:length(W.TF_names);
Wblock                = W;
TF_removed            = [];

% --------------------------------------------------
% ... or select data

% I. choose one block 

%it = 5;  
%gene_indices_choice = W.gene_sets{it};
%TF_indices_choice     = W.TF_sets{it};

% II. choose genes regulated by certain transcription factors

%[Wsub,gene_indices_choice,TF_indices_choice,Wsubpicture] = pick_subW_TF(W,{'crp'},1)

% III. Choose a few-input-module

%[n_TF,n_genes,TF_involved] = find_FIMS(W);
%dummy = find(n_TF==2);
%table(W.TF_names(dummy));
%ind = 3; % trpR tyrR
%genes = W.gene_names(find(W.data(:,dummy(ind))));
%[Wsub,gene_indices_choice,TF_indices_choice,Wsubpicture] = pick_subW_genes(W,genes,1);

% -----------------------------------------------------

if identifiable_flag,

  display(['Selecting identifiable subsystem']);

  [Wblock,genes_kept,TF_removed,genes_removed] = ...
      make_W_identifiable(W,gene_indices_choice,TF_indices_choice,0);
  
  gene_indices_choice = gene_indices_choice(genes_kept);
  
  Wblock.gene_names   = D.gene_names(gene_indices_choice);
  Wblock.operon_names = D.operon_names(gene_indices_choice);
  Wblock.operon_flags = D.operon_flags(gene_indices_choice);
  Wblock.operon_abbr  = D.operon_abbr(gene_indices_choice);
  
  [Wblock,gene_order,TF_order,Wpicture] = sort_W_matrix(Wblock);
  
  gene_indices_choice = gene_indices_choice(gene_order);
  DD = choose_from_D(D,gene_indices_choice);
  W = Wblock;
  D = DD;
else,
  display(['The system may be non-identifiable.']);  
end

display('Selecting TF and genes according to given list');

if isempty(gene_list), gene_list = W.gene_names; end
if isempty(TF_list),     TF_list = W.TF_names; end

[W,D] =  reorder_genes_W_D(gene_list,TF_list,W,D);

cd_data(filename,'analysis');
save(['identifiable_' filename],'D','W','TF_removed','gene_indices_choice'); 

display([' Writing file "identifiable_' filename '"']);
