function [n_TF,n_operons,TF_involved] = find_FIMS(Wblock)

% [n_TF,n_operons,TF_involved] = find_FIMS(Wblock)

for it=1:size(Wblock.data,2),
   n_operons(it) = sum(Wblock.data(:,it));
   TF_involved{it} =  find( sum( Wblock.data(find(Wblock.data(:,it)),:),1) ~=0);
   n_TF(it) = sum( sum( Wblock.data(find(Wblock.data(:,it)),:),1) ~=0);

end
