function [a,X_move,X_pred,m] = robust_scale_to_median(XX);

% [a,X_move,X_pred] = robust_scale_to_median(XX);

a = ones(1,size(XX,1));
for it2=1:10,
%m = nanmedian(repmat(a',1,size(XX,2)).*XX);
m = nanmean(repmat(a',1,size(XX,2)).*XX);
 for it = 1:size(XX,1), a(it) = median_match(m,XX(it,:)); end
end
 
X_pred =  (1./a)'*m;
X_move =  (repmat(a',1,size(XX,2)).*XX);
plot(XX(1,:)); hold on; plot(m,'r'); hold off