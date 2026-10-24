%-----------------------------------------------
% Application of NCA to artificial data example 
%-----------------------------------------------

clear

%------------------------------------------
% define true model

gene_names = {'X','Y','Z'}';
TF_names   = {'TF1','TF2'}';
A_true     = [1 0; 0.5 0.5; 0 -1]; 
B_true     = [sin(0:0.2:pi); cos(0:0.2:pi)];


%------------------------------------------
% nca input: network and artificial data 

W.data       = double(A_true~=0);
W.signs      = sign(A_true);
W.gene_names = gene_names;
W.TF_names   = TF_names;


%------------------------------------------
% nca input: artificial data 

X            = A_true * B_true;


%------------------------------------------
% run nca

[A,B,X_pred,info,A_list,B_list] = nca(X, W.data, W.signs,'graphics_flag',0,'verbose',1,'epsilon',10^-5);


%------------------------------------------
% display results as matrices

nca_display(X,A,B,A*B,W,1,12)


%------------------------------------------
% display network and curves

display_curves_and_connectivity(W,X,B);
