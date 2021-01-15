function [W2,Wp2,W_combi] = combine_networks(W,Wp,TF_rename_list,operon_rename_list,n_min,add_TF)

% [W1_new,W2_new,W_combined] = combine_networks(W1,W2,rename_list)
%
% Add entries from network W2 to W1.
% Compare the operons according to the field operon_names.
% Replace TF_names and operon_names in W2 according to TF_rename_list,operon_rename_list
% the pairs in these lists correspond to the names in W1 and W2, respectively.
% the names according to W1 are kept.
%
% match differing TF names according to rename_list
% (gene names are supposed to be from the same dictionary)

for it = 1:length(TF_rename_list),
 if find_TF(TF_rename_list{it}{2},Wp),
  Wp.TF_names{find_TF(TF_rename_list{it}{2},Wp)} = TF_rename_list{it}{1};
 end
end

for it = 1:length(operon_rename_list),
 if find_operon(operon_rename_list{it}{2},Wp);
  Wp.operon_names{find_operon(operon_rename_list{it}{2},Wp)} = operon_rename_list{it}{1};
 end
end

% --------------------------------------------------
% create a structure Wp1: same format as W, but with network entries 
% (data, signs, numbers) according to Wp

lg = label_names(W.operon_names,Wp.operon_names);
lt = label_names(W.TF_names,Wp.TF_names);

Wp1                            = W;
Wp1.data                       = zeros(size(W.data));
Wp1.signs                      = zeros(size(W.data));
Wp1.numbers                    = zeros(size(W.data));
Wp1.data(find(lg),find(lt))    = Wp.data(   lg(find(lg)),lt(find(lt)));
Wp1.signs(find(lg),find(lt))   = Wp.signs(  lg(find(lg)),lt(find(lt)));
Wp1.numbers(find(lg),find(lt)) = Wp.numbers(lg(find(lg)),lt(find(lt)));

% -------------------------------------------------
% also include transcription factors contained in Wp and not in W

if ~add_TF,

 W2  = W;
 Wp2 = Wp1;

else,

ltr = label_names(Wp.TF_names,W.TF_names);

if add_TF, new_TFs = find(ltr ==0);
else,  new_TFs = [];
end

relevant_new_TFs = new_TFs(find(sum(Wp.numbers(lg(find(lg)),new_TFs) >= n_min)));

Wp2.data  = double(Wp1.numbers>=n_min);
Wp2.signs = Wp1.signs;

Wp2.TF_names = [Wp1.TF_names; Wp.TF_names(relevant_new_TFs)];
Wp2.data(find(lg),size(Wp1.data,2)+1:size(Wp1.data,2)+length(relevant_new_TFs)) = ...
    double(Wp.data(lg(find(lg)),relevant_new_TFs)>=n_min);
 Wp2.signs(find(lg),size(Wp1.data,2)+1:size(Wp1.data,2)+length(relevant_new_TFs)) = ...
     Wp.signs(lg(find(lg)),relevant_new_TFs);
 Wp2.numbers = Wp1.numbers;
 Wp2.numbers(find(lg),size(Wp1.numbers,2)+1:size(Wp1.numbers,2)+length(relevant_new_TFs)) = ...
     Wp.numbers(lg(find(lg)),relevant_new_TFs);
if isfield(Wp1,'gene_names'),  Wp2.gene_names   = Wp1.gene_names; end
 Wp2.operon_names = Wp1.operon_names;
 
if isfield(Wp1,'operon_abbr'), Wp2.operon_abbr = Wp1.operon_abbr; end

 W2.TF_names = [W.TF_names;  Wp.TF_names(relevant_new_TFs)];
if isfield(Wp1,'gene_names'),   W2.gene_names = Wp1.gene_names; end
 W2.operon_names = Wp1.operon_names;
if isfield(Wp1,'operon_abbr'),  W2.operon_abbr = Wp1.operon_abbr; end
 W2.data = W.data;
 W2.data(:,size(Wp1.data,2)+1:size(Wp1.data,2)+length(relevant_new_TFs)) = 0;
 W2.signs = W.signs;
 W2.signs(:,size(Wp1.data,2)+1:size(Wp1.data,2)+length(relevant_new_TFs)) = 0;

end

% -----------------------------------------------------
% make combined matrix containing the information of both W and Wp

if nargout ==3,

W_combi.data  = double(W2.data + Wp2.data>0);
W_combi.signs = W2.signs;
W_combi.signs((W2.signs==0)) =  Wp2.signs(W2.signs==0);
W_combi.signs(find((W2.signs ~= Wp2.signs) .* (W2.signs .* Wp2.signs ~=0)) ) =  0;

W_combi.TF_names     = W2.TF_names;
if isfield(W2,'gene_names'),  W_combi.gene_names   = W2.gene_names; end
W_combi.operon_names = W2.operon_names;
if isfield(Wp1,'operon_abbr'), W_combi.operon_abbr  = W2.operon_abbr; end

end
