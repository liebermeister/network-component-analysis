function msr = hill_msr(kpar,x,y)

ymax = kpar(1);
km = kpar(2);
hc = kpar(3);

y_pred = hill(x,ymax,km,hc);
msr = mean((y-y_pred).^2);