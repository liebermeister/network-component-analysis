function [W, A, weight_quantity_type] = nca_network_sbtab_load(input, options)

% NCA_NETWORK_SBTAB_LOAD Read a transcription network (arrows and signs) from an SBtab file
%
% [W, A, weight_quantity_type] = nca_network_sbtab_load(input, options)
%
% Reads genes (Gene table), transcription factors (Regulator table) and the regulation
% arrows with their signs, by default from a Regulation table (SBtab 1.2), or from a
% column !ActivityCondition of the Gene table (e.g. '+ TF1 - TF2'; options.source =
% 'ActivityCondition'). Weights of the arrows (NCA matrix A) are read from a SparseMatrix
% table (see 'sbtab_regulation_network_load').
%
% input    SBtab document object or filename
% options  see 'sbtab_regulation_network_load' (fields source, target_table_type, ...)
%
% Output: network structure W (see 'help nca_networks') with fields
%   gene_names, TF_names, data (genes x TFs: 1 for each arrow),
%   signs (+1 activation, -1 inhibition, 0 unknown or dual), operon_names (= gene_names)
% A: weights (genes x TFs; [] if no weights table is given)
% weight_quantity_type: quantity type of the weights, e.g. 'regulation weight (log scale)'
%
% Requires the SBtab toolbox (https://github.com/liebermeister/sbtab-matlab).
%
% See also: nca_network_sbtab_save, sbtab_regulation_network_load

eval(default('options','struct'));
if isempty(options), options = struct; end

network = sbtab_regulation_network_load(input, options);

W.gene_names   = network.targets;
W.TF_names     = network.regulators;
W.data         = full(network.arrows);
W.signs        = full(network.signs);
W.operon_names = W.gene_names;

A = []; weight_quantity_type = '';
if isfield(network,'weights')
  A = full(network.weights);
  weight_quantity_type = network.weight_quantity_type;
end
