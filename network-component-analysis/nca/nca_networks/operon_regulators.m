function t = operon_regulators(operon,W)

% t = operon_regulators(operon,W)
% find the regulators of a operon

% t = regulators(operon,W)
%

t = W.TF_names( find(W.data(find_operon(operon,W),:)));
