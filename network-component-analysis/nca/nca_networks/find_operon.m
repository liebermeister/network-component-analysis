function index = find_operon(operon,D)

if isstr(operon),
index = label_names({operon},D.operon_names,'multiple');
index = index{1};
else,
index = label_names(operon,D.operon_names);
end

if isempty(index), fprintf('operon %s not found\n',operon); end