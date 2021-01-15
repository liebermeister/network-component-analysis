function gfp_std = local_time_shift(gfp_mat,od,std_od)

maxshift = 20;

for it =1:size(gfp_mat,1),

gfp = gfp_mat(it,:);

clear dum
for it2 = 1:2*maxshift,
  dum(it2,:) = circshift([std_od nan *ones(1,2*maxshift)],[0 it2]);
end
dum = dum - repmat([nan *ones(1,maxshift) od(it,:) nan *ones(1,maxshift)],2*maxshift,1);

[dummy,shifts] = min(abs(dum));
shifts = shifts - maxshift;
shifts(1:maxshift+4)= shifts(maxshift+5);
shifts(end-maxshift+1:end)= shifts(end-maxshift);
shifts =  0.5*conv2(shifts,1/maxshift*ones(1,maxshift),'same');

%plot(shifts)

gfp = fill_nan(gfp);
gfp_std(it,:) = interp1( (1:length(gfp) ) + shifts(maxshift+1:end-maxshift) , gfp, 1:length(gfp) ,'pchip'); 
gfp_std(it,end-ceil(shifts(maxshift+1:end-maxshift)):end) = gfp(end);

end

%figure(1); plot(gfp_std); hold on; plot(gfp,'r'); hold off
%figure(2); plot(std_od); hold on; plot(od,'r'); hold off
