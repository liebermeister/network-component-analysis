function y = hill(x,ymax,km,hc)

y = ymax ./ (1 + (x./km).^-hc);
y(x<0)=0;