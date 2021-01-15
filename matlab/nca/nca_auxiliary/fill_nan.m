function M= fill_nan(M)

for it = 1:size(M,1),
   dum = find(isfinite(M(it,:)));
   M(it,1:min(dum)-1)=M(it,min(dum));
   M(it,max(dum)+1:end)=M(it,max(dum));
end
