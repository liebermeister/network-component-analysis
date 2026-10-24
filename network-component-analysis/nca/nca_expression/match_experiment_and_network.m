function [W,D,W_picture,D_unmatched] = match_experiment_and_network(network_name,data_name,version_name);

% match_experiment_and_network(network_name,data_name,version_name);
% map data matrix to regulation matrix
% reorder promoters and transcription factors according to the given lists 

cd_data(data_name,'expression');
load(['data_' data_name ]);

cd_data(network_name,'network');
load(['network_' network_name]);

display(['Matching data "' data_name '" and network "' network_name '"']);

% --------------------------------------------

%%%%%% NOT NECESSARY FOR AMINO ACID DATA

if ~isfield(D,'operon_names'), 
  [D.operon_names,D.operon_flags,D.operon_abbr] = e_coli_gene2operon(D.gene_names);
end

[D_comb,W_comb,W_picture,D_unmatched] = combine_D_and_W(D,W);

% figure(3); D_display(D_comb);

% --------------------------------------------

W = W_comb;
D = D_comb;

filename = [data_name '__' network_name version_name];

cd_data(filename,'analysis');
save(['matched_' filename], 'W','D','W_picture','D_unmatched');
display([' Writing file "matched_' filename '"']);
