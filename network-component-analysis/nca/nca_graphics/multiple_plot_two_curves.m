function multiple_plot_two_curves(B1,B2,names,s1,s2,B3,s3);
% multiple_plot_two_curves(B1,B2,s1,s2,names,B3,s3);

if ~exist('s1','var'), s1 = 'b'; end
if ~exist('s2','var'), s2 = 'r--'; end
if ~exist('s3','var'), s3 = 'g'; end

[ni,nk] = subplot_n(size(B1,1));
for it = 1:size(B1,1),
   subplot(nk,ni,it); plot(B1(it,:),s1); hold on; plot(B2(it,:),s2); hold off
 if exist('names','var'), xlabel(names{it}); end
% noticks
%set(gca,'YScale','log');

end

if nargin>5
for it = 1:size(B1,1),
   subplot(nk,ni,it); hold on; plot(B3(it,:),s3); hold off
end

end
