function smooth_different(x,x_error)

x_start=X(3,:);
x=x_start;
x_error = x_error/(2*max(x_error));


for it = 1:1000
 x = [x(1) ...
      x(2:end-1).*(1-x_error(2:end-1)) + 0.5*x_error(2:end-1).* ( x(3:end) + x(1:end-2)),...
      x(end)];
plot(x_start,'r');
 hold on; plot(x,'g'); hold off
pause;
end
