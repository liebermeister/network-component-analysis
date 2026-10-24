function multiple_plot(M,names,error)

[ni,nk] = subplot_n(size(M,1));

for it =1 : size(M,1),
  subplot(ni,nk,it);
  plot(M(it,:),'r'); if exist('names','var'),title(names{it}); end
  axis tight
end

if exist('error','var'),
  if length(error)==1,
    error = error * ones(size(M));
  end
  for it =1 : size(M,1),
    subplot(ni,nk,it);
    hold on
    plot(M(it,:)+error(it,:),'b');
    plot(M(it,:)-error(it,:),'b');
    axis tight
  end
end
