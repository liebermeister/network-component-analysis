function [p_value,sign_change,A_median_1,A_median_2] = nca_significant_A_differences(A_list1,A_list2,W)

% determine where two groups of A matrices (e.g., resulting from two sift experiments) show significant differences
%
% attention: all matrices need to have the same format (= correspond to the same W matrix)
%
% p values are computed by 'ranksum', i.e., a rank test for different median

for it = 1:length(A_list1),
  A_tensor1(:,:,it) = A_list1{it};
end

for it = 1:length(A_list1),
  A_tensor2(:,:,it) = A_list2{it};
end

A_median_1  = median(A_tensor1,3);
A_median_2  = median(A_tensor2,3);

ind_arrows = find(W.data);
[ind_gene,ind_TF] = ind2sub(size(W.data),ind_arrows);

p_value = nan * A_median_1;
for it = 1:length(ind_arrows),
  p_value(ind_gene(it),ind_TF(it)) = ranksum(squeeze(A_tensor1(ind_gene(it),ind_TF(it),:)),squeeze(A_tensor2(ind_gene(it),ind_TF(it),:)));
end

sign_change = [sign(A_median_1) ~= sign(A_median_2)];
