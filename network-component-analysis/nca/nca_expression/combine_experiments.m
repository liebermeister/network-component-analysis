function D3 = combine_experiments(D1,D2,name1,name2)

% D3 = combine_experiments(D1,D2)
%
% infos in D2 override those in D1
% gene names must not appear multiply within D1 or D2
% gene names must match between D1 and D2

gene_names = [setdiff(D1.gene_names,D2.gene_names); ...
              intersect(D1.gene_names,D2.gene_names); ...
              setdiff(D2.gene_names,D1.gene_names)];

l1 = label_names(D1.gene_names,gene_names);
l2 = label_names(D2.gene_names,gene_names);

D3.gene_names = gene_names;

if isfield(D1,'gene_annotations') &  isfield(D2,'gene_annotations'),
 D3.gene_annotations(l1) = D1.gene_annotations;
 D3.gene_annotations(l2) = D2.gene_annotations;
end

lt1 = size(D1.GFP_shifted,2);
lt2 = size(D2.GFP_shifted,2);

if ~isfield(D1,'experiments'), D1.experiments = {1:lt1}; D1.experiment_names{1}=name1; end
if ~isfield(D2,'experiments'), D2.experiments = {(1:lt2)+lt1}; D2.experiment_names{1}=name2; end

D3.GFP_shifted = nan*ones(length(gene_names), lt1+lt2);
D3.GFP_shifted(l1,1:lt1) = D1.GFP_shifted;
D3.GFP_shifted(l2,lt1+1:end) = D2.GFP_shifted;

D3.OD_shifted = nan*ones(length(gene_names), lt1+lt2);
D3.OD_shifted(l1,1:lt1) = D1.OD_shifted;
D3.OD_shifted(l2,lt1+1:end) = D2.OD_shifted;

D3.GFP_p_OD = nan*ones(length(gene_names), lt1+lt2);
D3.GFP_p_OD(l1,1:lt1) = D1.GFP_p_OD;
D3.GFP_p_OD(l2,lt1+1:end) = D2.GFP_p_OD;

D3.GFPder_p_OD = nan*ones(length(gene_names), lt1+lt2);
D3.GFPder_p_OD(l1,1:lt1) = D1.GFPder_p_OD;
D3.GFPder_p_OD(l2,lt1+1:end) = D2.GFPder_p_OD;

D3.growth_rate = nan*ones(length(gene_names), lt1+lt2);
D3.growth_rate(l1,1:lt1) = D1.growth_rate;
D3.growth_rate(l2,lt1+1:end) = D2.growth_rate;

D3.GFPder_p_OD_Std_err = nan*ones(length(gene_names), lt1+lt2);
D3.GFPder_p_OD_Std_err(l1,1:lt1) = D1.GFPder_p_OD_Std_err;
D3.GFPder_p_OD_Std_err(l2,lt1+1:end) = D2.GFPder_p_OD_Std_err;

D3.experiment_names  = [D1.experiment_names, D2.experiment_names];
D3.experiments = D1.experiments;
for it=1:length(D2.experiments),
  D3.experiments = [D3.experiments, {D2.experiments{it}+size(D1.GFP_shifted,2)}];
end;

