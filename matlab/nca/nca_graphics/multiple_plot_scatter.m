function multiple_plot_scatter(B1,B2);

[ni,nk] = subplot_n(size(B1,1));
for it = 1:size(B1,1),
   subplot(ni,nk,it); plot(B1(it,:), B2(it,:),'.');
end
