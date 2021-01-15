function [Wsub,operon_indices,TF_indices,Wpicture] = pick_subW_TF(W,TFlist,sorting,omit)

% [Wsub,operon_indices,TF_indices,Wpicture] = pick_subW_TF(W,TFlist,sorting,omit)

if ~exist('sorting','var'), sorting = 1; end
if ~exist('omit','var'),    omit    = 0; end
if ~exist('strict','var'),  strict  = 0; end

TFindices = find_TF(TFlist,W);
TFindices = TFindices(find(TFindices)); 

if strcmp(omit,'direct'),

  TF_indices     = TFindices;
  operon_indices = find(sum(W.data(:,TF_indices),2));

else,

  operon_indices = find(sum(W.data(:,TFindices),2));
  TF_indices     = find(sum(W.data(operon_indices,:),1));
  
  if strict,
    TF_indices     = TFindices;
    operon_indices = find((sum(W.data(:,TFindices),2)) .* ( sum(W.data,2) <= sum(W.data(:,TFindices),2)));
  end
  
  if length(TF_indices)==1,
    if omit,
      TF_indices = setdiff(TF_indices,TFindices); 
      fprintf('Omitting common transcription factor %s\n',W.TF_names{TFindices});
    end
end

end

Wsub.data         = W.data(operon_indices,TF_indices);
Wsub.signs        = W.signs(operon_indices,TF_indices);
Wsub.operon_names = W.operon_names(operon_indices);
if isfield(W,'gene_names'),  
  Wsub.gene_names=W.gene_names(operon_indices); 
end
Wsub.TF_names     = W.TF_names(TF_indices);

Wpicture = Wsub.data;
if sorting, [Wsub,Wpicture] = sort_W_matrix(Wsub); end
