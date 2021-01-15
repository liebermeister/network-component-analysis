% Combine all data that form an identifiable network
% with all data that do not at all match the regDB network
%
% prepare the data structure D and according network structures describing the
% known and putative binding sites

function care_for_unmatched_genes(filename,TF_list,gene_list,A,sigma_flag)

if sigma_flag == 0, 
  TF_list = TF_list(find(strncmp('Sigma',TF_list,5)==0));
end


% ---------------------------------------
% get putative network for all unmatched genes: Wp_unmatched 

cd_data(filename,'analysis');
load(['matched_' filename]); 

lTF = label_names(TF_list,W.TF_names);
l_gene = label_names(gene_list,W.gene_names);

W = choose_from_W(W,l_gene,lTF);
D = choose_from_D(D,l_gene);

%get D, D_unmatched, W 

cd_data('','network');
load network_putative;

[D2,Wp_unmatched] = combine_D_and_W(D_unmatched,Wp);

% ---------------------------------------
% get matched identifiable network: D and W 

cd_data(filename,'analysis');
load([ 'identifiable_' filename ]);
load([ 'final_' filename ]);

lTF = label_names(TF_list,W.TF_names);
l_gene = label_names(gene_list,W.gene_names);

W = choose_from_W(W,l_gene,lTF);
D = choose_from_D(D,l_gene);

% ---------------------------------------

W1.TF_names     = W.TF_names; 
W1.gene_names   = [W.gene_names; Wp_unmatched.gene_names];
W1.operon_names = [W.operon_names; Wp_unmatched.operon_names];
W1.operon_abbr  = [W.operon_abbr; Wp_unmatched.operon_abbr];
W1.operon_flags = [W.operon_flags; Wp_unmatched.operon_flags];
W1.data         = [W.data; zeros(length(Wp_unmatched.operon_names),size(W.data,2))];
W1.signs        = [W.signs; zeros(length(Wp_unmatched.operon_names),size(W.data,2))];

%rename_list = {{'ompR envZ','ompR'},{'fis-dusB','fis'},{'nagBACD','nagC'}};
rename_list = {};
[dummy1,dummy2,W_combi] = combine_networks(W1,Wp,rename_list,{},1,0);

D_combi.operon_names        = [D.operon_names; D2.operon_names];
D_combi.gene_names          = [D.gene_names; D2.gene_names];
D_combi.GFP_p_OD            = [D.GFP_p_OD; D2.GFP_p_OD];
D_combi.GFPder_p_OD         = [D.GFPder_p_OD; D2.GFPder_p_OD];
D_combi.GFP_shifted         = [D.GFP_shifted; D2.GFP_shifted];
D_combi.OD_shifted          = [D.OD_shifted; D2.OD_shifted];
D_combi.GFPder_p_OD_Std_err = [D.GFPder_p_OD_Std_err; D2.GFPder_p_OD_Std_err];
D_combi.growth_rate         = [D.growth_rate; D2.growth_rate];
D_combi.experiments         = D.experiments;
D_combi.experiment_names    = D.experiment_names;
D_combi.time_point_names    = D.time_point_names;

blocks = D.experiments;
X_combi      = [];
X_combi_err  = [];
for it = 1:length(blocks),
 X_combi= [X_combi fill_nan(D_combi.GFPder_p_OD(:,blocks{it}))];
 this_Xerr = D_combi.GFPder_p_OD_Std_err(:,blocks{it});
 this_Xerr(isnan(D_combi.GFPder_p_OD(:,blocks{it}))) = 300;
 X_combi_err = [X_combi_err  this_Xerr];
end
  
X_combiplus = X_combi + X_combi_err;
X_combi = log(X_combi+100)-log(100);
X_combiplus = log(X_combiplus+100)-log(100);
X_combi_error = X_combiplus-X_combi;

%figure(1); display_curves_and_connectivity(W,X);
%figure(2); display_curves_and_connectivity(W1,X_combi);
%figure(3); display_curves_and_connectivity(W_combi,X_combi);

A_combi = [A; zeros(length(D_unmatched.operon_names),size(A,2))];

save([filename '_plus_unmatched_putative'],'W1','W_combi','D_combi','X_combi','X_combi_error','A_combi');
