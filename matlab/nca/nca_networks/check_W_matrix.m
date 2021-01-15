function check_W_matrix(W)

% check_W_matrix(W)
% little test for network data structure: find lac operon and targets of crp

lacz_ind            = find(strcmp(W.operon_names,'lacZYA'));
W.operon_names{lacz_ind}
table(W.TF_names(find(W.data(lacz_ind,:))))
crp_ind = find(strcmp(W.TF_names,'lacI'));
W.TF_names{crp_ind}
table(W.operon_names(find(W.data(:,crp_ind))));

