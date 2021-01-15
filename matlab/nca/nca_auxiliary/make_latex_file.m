function text = make_latex_file(insert_text,file_name)

starttext = sprintf('\\documentclass[a4paper]{article}\n\\usepackage{setspace}\n\\usepackage{epsfig}\n\\usepackage{graphicx}\n\\begin{document}\n');

endtext = sprintf('\\end{document}');

text = [starttext, insert_text, endtext];

if exist('file_name','var'),
  fid = fopen([file_name '.tex'],'wt');
  fprintf(fid,'%s',text);
  fclose(fid);

  eval(['! sh -c "/usr/local/package/bin/latex ' file_name '.tex"']);
  eval(['! sh -c "/usr/local/package/bin/dvips ' file_name '.dvi"']);
end

