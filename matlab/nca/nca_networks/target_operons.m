function t = target_operons(TF,W)

% t = targets(TF,W)
% find the target operons of a TF 
% examples: 
%  target_operons('argR',W)
%  target_operons({'argR','crp'},W)

if iscell(TF),
    t = {};
    for it =1:length(TF),
 t = [t; W.operon_names( find(W.data(:,find_TF(TF{it},W))))];
end
    t = unique(t);
else
  t = W.operon_names( find(W.data(:,find_TF(TF,W))));
end
