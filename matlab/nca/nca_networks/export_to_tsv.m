function T = export_to_tsv(W,fields,headers,filename,A,weightnames)

% export_to_tsv(W,fields,headers,filename)
%
% filename is optional

eval(default('headers','{''gene_names'',''TF_names''}'));
eval(default('fields','{''gene_names'',''TF_names''}'));
eval(default('filename','[]'));
eval(default('weightnames','''Weight'''));

[row,col] = ind2sub(size(W.data),find(W.data));

n1 = getfield(W,fields{1});
n2 = getfield(W,fields{2});

signs = W.signs(find(W.data));

signsstr = {};
signsstr(signs==1,1) = {'+1'};
signsstr(signs==0,1) = {'unknown'};
signsstr(signs==-1,1) = {'-1'};

T = [[headers {'Sign'}]; [n1(row) n2(col) signsstr]];

if exist('A','var'),
  if ~iscell(A),           A           = {A};           end
  if ~iscell(weightnames), weightnames = {weightnames}; end
  weightnames = column(weightnames)';
  for it = 1:length(A),
    T = [ T, [weightnames(it); num2cellstr(A{it}(find(W.data)))]];
  end
end

if ~isempty(filename), mytable(T,0,filename); end