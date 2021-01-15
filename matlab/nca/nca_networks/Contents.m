%Data structure for gene regulatory networks used in NCA
%
%Wolfram Liebermeister, 2006
%
%A network with nt transcription factors and np promoters is represented by a data structure
%
%DATA STRUCTURE:
%
%Mandatory fields:
% gene_names:          list of strings
% TF_names:            list of strings
% data  (np x nt):     0: no interaction, 1:interaction
% signs (np x nt):     0: no/dual interaction, 1:activation, -1:inhibition
%
%Optional fields:
% operon_names        list of strings: operon names
% operon_abbr         short names of operons
% operon_flags        (1: gene is first gene in operon;  0: not first gene in operon;  -1: operon not found)
% gene_sets
% TF_sets
% numbers
%
%FUNCTIONS:
%
%Find things in network W
%
%  find_operon('argE',W)                find index of an operon (same for 'find_gene')
%  find_operon({'argE','argI'},W)       same, for list of operons (same for 'find_gene')
%  find_TF('crp',W)                     find index of a transcription factor
%  find_TF({'crp','lrp'},W)             same, for list of transcription factors
%  find_FIMS                            find small subnetworks
%  find_FIMS(Wblock)                    modules (one TF + its targets+their regulators)
%  regulators(gene,W)                   find the regulators of a gene
%  target_genes(TF,W)                   find the target genes of a TF 
%  target_operons(TF,W)                 find the target operons of a TF 
%
%Manipulate network
%
%  sort_W_matrix                        sort network structure -> block structure
%  pick_subW_TF(W,'crp',0)              Choose subnetwork
%  pick_subW_operons(W,'argE')          Choose subnetwork
%  combine_networks                     merge two TF networks
%  rearrange_W                          reorder network/ pick subnetwork
%  W_add_general                        add a general TF
%  randomise_edges                      create a random network with same degrees
%
%Display network
%
%  W_display(Wsub);                     display network connectivities by lines
%  W_display(Wsub,'N');                 display network connectivities as a matrix
%
%Display network together with expression data
%
%  display_curves_and_connectivity      display network connectivities, add time courses
%  display_D_block                      display network and experimental data
