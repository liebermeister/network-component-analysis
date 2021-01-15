% Scripts to preprocess expression data for NCA
%
% How to turn matrices called 'od' and 'gfp' into a data structure 'D'
% This structure contains, among other things, smooth curves for 
% GFP/OD and (dGFP/dt)/OD
%
%1. put data into a structure D
%  o each matrix row of gfp and od must correspond to one gene
%  o matrices are supposed to be already preprocessed by alon's matlab script
%  o gene_names is a (column) list of strings
%
% D.GFP        = gfp;
% D.OD         = od;
% D.gene_names = gene_names;
%
%2. run the whole preprocessing -> updated structure D
%
% [D,info] = preprocess_gfp_data(D,'graphics_flag',0,'remove_genes_flag',1,'use_time_points',21:size(D.OD,2));
%
%3. If necessary, average over multiple occurrences of genes 
% D = gfp_data_unify_genes(D);
%   (this can also be done automatically within 'preprocess_gfp_data')
%
% For more information, type 'help process_gfp_data'
%II. Expression Data Functions
%
%Experimental data and transcriptional regulation networks are
%represented by matlab structures.
%
%Mandatory fields for experimental data structure D
% operon_names  list of (m) operon names
% GFPder_p_OD   m x k matrix containing (d GFP/dt)/OD values for k time points 
%
%Mandatory fields for transcriptional networks structure W
% operon_names  list of (m) operon names
% TF_names      list of (n) transcription factor names
% data          m x n matrix containing 1 if regulation is known, 0 otherwise
% signs         similar, with 1: positive control. 2: negative control. 3: dual control
%
%Manipulate Data Matrices
%
%  clean_gfp_data                 estimate gfp'/od etc
%  check_control_genes            plot timecourses of repeats
%  combine_D_and_W                make a data structure and a network structure comparable
%  combine_experiments            combine two data structures
