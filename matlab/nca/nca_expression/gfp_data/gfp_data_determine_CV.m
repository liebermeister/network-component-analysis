function mean_cv = gfp_data_determine_CV(D);

gene_names = unique(D.gene_names);

l = label_names(gene_names,D.gene_names,'multiple');

cv=[];
for it = 1:length(l),
  indices = l{it};
  if length(indices)>1,
%      length(indices)
%       D.gene_names(l{it}(1))
      gfp_mean = mean(D.GFPder_p_OD(indices,:));
      gfp_std = std(D.GFPder_p_OD(indices,:));
      cv = [cv; gfp_std ./ (gfp_mean+100)];
end
end

%plot(cv'); hold on;
%plot(mean(cv),'*'); hold off
mean_cv = mean(cv);
