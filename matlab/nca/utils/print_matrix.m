function b=print_matrix(matrix,rownames,colnames,number_flag)

% b=print_matrix(matrix,rownames,colnames,number_flag)

if ~exist('number_flag','var'), number_flag = 0; end

if ~exist('rownames','var'), rownames=cell(size(matrix,1),1); end
if ~exist('colnames','var'), colnames={}; end
  
  if isempty(colnames),
 b= [ char(rownames) repmat(' ',size(matrix,1),1) num2str(matrix)];

 if number_flag, b=[ num2str((1:size(matrix,1))') repmat(' ',size(matrix,1),1) b]; end 
  else,
  
if size(rownames,1) ==1,  rownames=rownames'; end
if size(colnames,1) ==1,  colnames=colnames'; end

 a=num2cell(matrix);
 b=cell(size(a,1)+1,size(a,2)+1);
 b(2:end,2:end)=a;
 b(1,2:end)=colnames';
 b(2:end,1)=rownames;
end