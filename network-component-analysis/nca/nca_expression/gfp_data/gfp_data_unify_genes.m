function DD = gfp_data_unify_genes(D);

% DD = gfp_data_unify_genes(D);
%
% identify genes that occur more than one time in the data set
% shift + average the respective time courses

DD.gene_names = unique(D.gene_names);

l = label_names(DD.gene_names,D.gene_names,'multiple');

for it = 1:length(l),
  indices = l{it};
  DD.original_indices{it} = D.original_indices(indices);
  DD.gene_annotations{it} = D.gene_annotations(indices(1));
  mediancurve = nanmedian(D.GFPder_p_OD(indices,:));
  clear GFPder_p_OD_shifted

  if length(indices)>1,
   for it2 = 1:length(indices);
    shift = transform_OD_timecourse(D.GFPder_p_OD(indices(it2),:)',mediancurve',[],10);
    this_GFP_shifted(it2,:) = shift_matrix(D.GFP_shifted(indices(it2),:),shift);
    this_OD_shifted(it2,:)  = shift_matrix(D.OD_shifted(indices(it2),:),shift);
    this_GFP_p_OD(it2,:)    = shift_matrix(D.GFP_p_OD(indices(it2),:),shift);
    this_GFPder_p_OD(it2,:) = shift_matrix(D.GFPder_p_OD(indices(it2),:),shift);
    this_growth_rate(it2,:) = shift_matrix(D.growth_rate(indices(it2),:),shift);
    if isfield(D,'GFPder_p_OD_Std_err'),
      this_GFPder_p_OD_Std_err(it2,:) = shift_matrix(D.GFPder_p_OD_Std_err(indices(it2),:),shift);
    end
   end
   DD.GFP_shifted(it,:) = nanmean(this_GFP_shifted);
   DD.OD_shifted(it,:)  = nanmean(this_OD_shifted);
   DD.GFP_p_OD(it,:)    = nanmean(this_GFP_p_OD);
   DD.GFPder_p_OD(it,:) = nanmean(this_GFPder_p_OD);
   DD.growth_rate(it,:) = nanmean(this_growth_rate);
   if isfield(D,'GFPder_p_OD_Std_err'),
     DD.GFPder_p_OD_Std_err(it,:) = nanmean(this_GFPder_p_OD_Std_err);
   end
   DD.shifts(it) = mean(D.shifts(indices));
  else,
   DD.GFP_shifted(it,:) = D.GFP_shifted(indices,:);
   DD.OD_shifted(it,:) = D.OD_shifted(indices,:);
   DD.GFP_p_OD(it,:)    = D.GFP_p_OD(indices,:);
   DD.GFPder_p_OD(it,:) = D.GFPder_p_OD(indices,:);
   DD.growth_rate(it,:) = D.growth_rate(indices,:);
   if isfield(D,'GFPder_p_OD_Std_err'),
     DD.GFPder_p_OD_Std_err(it,:) = D.GFPder_p_OD_Std_err(indices,:);
   end
   DD.shifts(it,:)      = D.shifts(indices);
  end
end
