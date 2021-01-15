function index = find_TF(TF,D)

% index = find_TF(TF,D)

if isstr(TF),
 index = label_names({TF},D.TF_names,'multiple');
 index = index{1};
else,
 index = label_names(TF,D.TF_names);
end

if isempty(index), fprintf('TF %s not found\n',TF); end