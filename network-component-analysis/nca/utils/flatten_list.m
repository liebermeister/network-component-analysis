function [items,occurrences,numbers] = flatten_list(list)

% takes a list of string lists and extracts all strings present
% occurrences is a list of index vectors describing where each item was found
% number is a list of indices corresponding to the original list, saying which
% items the list elements contain

items = [];

for it=1:length(list),
  if size(list{it},1)==1, list{it}=list{it}'; end
  items = [items; list{it}];
end

items = unique(items);


numbers = [];
occurrences = cell(size(items));
for it=1:length(list),
  numbers{it} = label_names(list{it},items,'single');
  for it2=1:length(numbers{it}),
    occurrences{numbers{it}(it2)} =  [  occurrences{numbers{it}(it2)}; it]; end
end

numbers = numbers';
