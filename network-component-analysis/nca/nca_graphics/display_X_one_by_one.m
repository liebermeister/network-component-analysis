if exist('X','var'),

for ind = 1:size(X,1),
 plot(X(ind,:))
 hold on; plot(X(ind,:)+X_error(ind,:),'r');
 hold on; plot(X(ind,:)-X_error(ind,:),'r'); hold off
pause
end

else,

for ind = 1:size(D.GFPder_p_OD,1),
 plot(D.GFPder_p_OD(ind,:))
 hold on; plot(D.GFPder_p_OD(ind,:)+D.GFPder_p_OD_Std_err(ind,:),'r');
 hold on; plot(max(0,D.GFPder_p_OD(ind,:)-D.GFPder_p_OD_Std_err(ind,:)),'r'); hold off
pause
end
end
