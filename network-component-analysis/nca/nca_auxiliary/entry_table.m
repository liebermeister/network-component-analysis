function entry_table(W, foundflag)


if exist('foundflag','var'),  
  found = find(W.operon_flags~=-1);
table([W.gene_names(found) cellstr(num2str(W.operon_flags(found))) W.operon_abbr(found) W.operon_names(found)])
else
 mytable([W.gene_names cellstr(num2str(W.operon_flags)) W.operon_abbr W.operon_names])
end