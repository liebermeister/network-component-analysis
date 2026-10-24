function changes_table(changes,W);

if length(changes.sign_changes),  
  sgns = diag(W.signs(changes.sign_changes(:,2),changes.sign_changes(:,1)));
  fprintf('Sign restrictions relaxed:\n');
 mytable( [ cellstr(num2str(sgns)) W.TF_names(changes.sign_changes(:,1)) W.operon_names(changes.sign_changes(:,2)) ]); 
end
if length(changes.new), 
  fprintf('Connections added:\n');
 mytable( [ W.TF_names(changes.new(:,1)) W.operon_names(changes.new(:,2)) ]);
end
if length(changes.removed), 
  fprintf('Connections removed:\n');
 mytable( [ W.TF_names(changes.removed(:,1)) W.operon_names(changes.removed(:,2)) ]);
end


if length(changes.new), 
  fprintf('Connections added:\n');
  names_string = repmat({''},size(W.data,2),1);
  for it =1:size(changes.new,1),
    this_TF = changes.new(it,1);
    if changes.new_sign(it)==1, sign_str = '+'; end
    if changes.new_sign(it)==-1, sign_str = '-'; end
    names_string{this_TF} = [ names_string{this_TF} ' ' sign_str W.operon_names{changes.new(it,2)}];
  end
  dummy = [ W.TF_names repmat({'&'},size(W.data,2),1) names_string  repmat({'\\'},size(W.data,2),1)];
 mytable(dummy(unique(changes.new(:,1)),:),0)
end
