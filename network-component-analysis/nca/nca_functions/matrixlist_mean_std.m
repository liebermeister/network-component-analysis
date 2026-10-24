function [Mmean,Mstd] = matrixlist_mean_std(Mlist,bfilt_flag,median_flag)

eval(default('bfilt_flag','0','median_flag','''median'''));

for it = 1:length(Mlist),
  if bfilt_flag,
    X_tensor(:,:,it)      = my_bfilt(Mlist{it});
  else,
    X_tensor(:,:,it)      = Mlist{it};
  end
end

switch median_flag,
  case 'median',
    for it1 = 1:size(X_tensor,1),
      for it2 = 1:size(X_tensor,2),
        c = column(squeeze(X_tensor(it1,it2,:)));
        Mmean(it1,it2) = nanmedian(c);
        Mstd(it1,it2)  = 0.5 *[quantile(c,0.75) - quantile(c,0.25)];
      end
    end
  otherwise,
    Mmean = squeeze(nanmean(X_tensor,3));
    Mstd  = squeeze(std(X_tensor,[],3));  
end
