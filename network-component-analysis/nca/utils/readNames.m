% lines = readNames(filename, split, nlines)
%
% split = 1 -> return cell structures containing tab-separated items
% split = 0 (default) return whole line as string

function names = readNames(filename,split,nlines)

file = fopen(filename);

if ~exist('split'), split=0; end
if exist('nlines'), names=cell(1,nlines); end

k=0; % line number

while(~feof(file))
  k = k+1;
  s = fgetl(file );
  if isstr(s)
   names{k} = deblank(s);
   if split
    nn = names{k};
    dum = [findstr(nn,char(9)) length(nn)+1 ];
    l={};
    if length(nn)>0
     l{1} = nn(1:dum(1)-1);
     for i=1:length(dum)-1
       l{i+1}=deblank(nn(dum(i)+1:dum(i+1)-1));
     end
    end
    names{k}=l;
   end
  end
end

fclose(file);

if length(names)==0, fprintf('readNames.m: Empty file %s\n', filename); end

names=names';
