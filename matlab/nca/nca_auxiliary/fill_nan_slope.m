function M= fill_nan_slope(M,xmin)

for it = 1:size(M,1),
   dum = find(isfinite(M(it,:)));
   M(it,1:min(dum)-1)=M(it,min(dum));
   slope = (M(it,dum(end))-M(it,dum(end-5))) /( dum(end) - dum(end-5) );
   M(it,max(dum)+1:end) =  max(xmin,M(it,dum(end)) + slope * (1:(size(M,2)-max(dum))));
end
