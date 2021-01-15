function X_display(X, operon_list,D,W,ni,color,heights)

if ~exist('ni','var'), ni = []; end
if isempty(ni), ni = length(operon_list); end
if ~exist('color','var'), color = [ 0 0 1]; end

for it = 1:length(D.experiments),
  starting_points(it)= D.experiments{it}(1);
end

for it = 1:length(operon_list),
  subplot(ni,1,it);
  index =   find_operon(operon_list(it),W);
  mmax = max(X(index,:));
  if exist('heights','var'), mmax = heights(it); end
  mmin = min(X(index,:));
  hold on;
  plot(X(index,:)','color',color);  
  for it2 = 1:length(D.experiments),
    line( starting_points(it2) *[1 1], [mmin,mmax],'color',[0 0 0]);
  axis([0  size(X,2) mmin  mmax]); 
%  axis tight;
  end
  hold off;
  ylabel(operon_list(it)); set(gca,'XTick',[],'YTick',[]);
end
set(gca,'XTick',starting_points','XTickLabel',D.experiments);
