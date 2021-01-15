function display_multiple_curves(data)

% display_multiple_curves(data)
% plot quantiles of many curves

subplot(2,1,1); 

plot(data','c'); hold on ;
mediancurve = nanmedian(data);
lower = nanquantile(data,0.05);
upper = nanquantile(data,0.95);
%lowerq = nanquantile(data,0.25);
%upperq = nanquantile(data,0.75);
plot(mediancurve,'k'); plot(mediancurve,'k.');
%plot(lowerq,'m');
%plot(upperq,'m');
plot(lower,'k'); plot(lower,'k.'); 
plot(upper,'k'); plot(upper,'k.'); hold off ;

subplot(2,1,2); im(data);
