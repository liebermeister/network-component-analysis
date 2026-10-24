function plot_true_versus_fit(X,X_pred,W,indices,siz)

% plot_true_versus_fit(X,X_pred,W,indices,siz)

if ~exist('indices','var'),
  indices = 1:size(X,1);
end

method = 'hill';

clear ymax km hc X_fit;

if exist('siz','var'),
  ni=siz(1);
  nk=siz(2);
  else,
[ni,nk] = subplot_n(length(indices));
end

for it = 1:length(indices),
  ind_operon = indices(it);
  subplot(ni,nk,it); 
  
  switch method,
    case 'interpolate',
      plot(X_pred(ind_operon,:),X(ind_operon,:),'.','color',[0.7,0.7,0.85]);
      [d2,order] = sort(X_pred(ind_operon,:)+10^-8*rand(1,size(X,2)));
      d1 = d2(1) + (d2(end)-d2(1))* (0:0.001:1);
      d3 = interp1(d2, X(ind_operon,order),d1,'cubic');
      hold on; 
      d4 = conv2(d3,1/200*ones(1,200),'same');
      plot(d1(100:end-100),d4(100:end-100),'r','linewidth',3);
      hold off; axis equal;  title(W.operon_names{ind_operon});
      noticks;  
    case 'hill',
      x = X_pred(ind_operon,:);
      y = X(ind_operon,:);
  % FOR UNTRANSFORMED DATA:
%  y = exp(y + log(100))-100;
      ymax_start = 1.5*max(y); km_start = max(x); hc_start = 1;
      kpar = [ymax_start km_start hc_start];
      kparopt = fminsearch(@hill_msr,kpar,optimset('MaxFunEvals',1000),x,y);
      ymax(it) = kparopt(1); km(it) = kparopt(2); hc(it) = kparopt(3);
      y_fit = hill(x,ymax(it),km(it),hc(it));
      X_fit(it,:)= y_fit;
      plot(x,y,'.','color',[0.7,0.7,0.85]); hold on 
      plot(sort(x),sort(y_fit),'r.'); hold off
  end
  noticks
  axis tight
  title(W.operon_names{ind_operon});
  drawnow
end
