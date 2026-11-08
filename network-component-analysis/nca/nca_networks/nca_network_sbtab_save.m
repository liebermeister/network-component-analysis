function sbtab_document = nca_network_sbtab_save(W, filename, options, A)

% NCA_NETWORK_SBTAB_SAVE Write a transcription network (arrows and signs) to an SBtab file
%
% sbtab_document = nca_network_sbtab_save(W, filename, options, A)
%
% Writes the genes (Gene table), the transcription factors (Regulator table, with !Type
% 'transcriptional regulator') and the regulation arrows with their signs, by default in a
% Regulation table (SBtab 1.2), or in a column !ActivityCondition of the Gene table
% (options.format = 'ActivityCondition'). If A (NCA weights, genes x TFs) is given, the
% weights are written to a SparseMatrix table; their quantity type is given by
% options.weight_quantity_type ('regulation weight (log scale)' (default) or
% 'regulation weight (linear scale)', depending on the scale of the expression data).
%
% W         network structure (see 'help nca_networks'; fields gene_names, TF_names,
%           data, signs)
% filename  output file (.tsv); if empty, no file is written
% options   see 'sbtab_regulation_network_save' (fields format, document_name, ...)
%
% Requires the SBtab toolbox (https://github.com/liebermeister/sbtab-matlab).
%
% See also: nca_network_sbtab_load, sbtab_regulation_network_save

eval(default('filename','[]','options','struct','A','[]'));
if isempty(options), options = struct; end

network.targets    = W.gene_names;
network.regulators = W.TF_names;
network.arrows     = W.data;
network.signs      = W.signs;
if ~isempty(A), network.weights = A; end

sbtab_document = sbtab_regulation_network_save(network, filename, options);
