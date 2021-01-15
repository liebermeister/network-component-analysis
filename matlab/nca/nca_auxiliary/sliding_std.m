function [Mmean,Mstd] = sliding_std(M,experiments)

window = 5;

for it=1:size(M,1),
  for it3 = 1:length(experiments),
    for it2 = experiments{it3}
      Mmean(it,it2) = mean( M(it,max(experiments{it3}(1),it2-window):min(experiments{it3}(end),it2+window )));
    end
  end
end

for it=1:size(M,1),
  for it3 = 1:length(experiments),
    for it2 = experiments{it3}
     Mstd(it,it2) = std( M(it,max(experiments{it3}(1),it2-window):min(experiments{it3}(end),it2+window )) -  Mmean(it,max(experiments{it3}(1),it2-window):min(experiments{it3}(end),it2+window )));
    end
  end
end