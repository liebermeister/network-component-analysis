function [y_derivative,y_smoothed] = smooth_derivative(y,n_window,method,x)

% [y_derivative,y_smoothed] = smooth_derivative(y,n_window,method,x)
%
% Smoothen time series and compute a robust time derivative
% The function is computed by regression (linear, square, or cubic) 
% within a sliding window of size 1+2*n_window around the respective data point.
%
% y            row vector of function values. if y is a matrix, each row is treated separately.
% n_window     half width (# data points) of sliding window (default 5)
% method       'linear','square', 'cubic'
% x            (optional) row vector of time points corresponding to function values in y
%              the values must be increasing
%              the default for x is equally spaced values 1..n 
%
% y_derivative estimated time derivative
% y_smoothed   smoothed version of y
%
%Test case 1:
% x = [0:0.1:10];
% y_true = sin(x) ;
% y = y_true + 0.1*randn(size(y_true));
% y_der_naive = [y(:,2)-y(:,1), 0.5*[y(:,3:end)-y(:,1:end-2)], y(:,end)-y(:,end-1)]/(x(1,2)-x(1,1));
% [y_derivative,y_smoothed] = smooth_derivative(y,20,'cubic',x);
% plot(x(:),y_true(:),'.'); hold on; 
% plot(x(:),y(:),'c.');
% plot(x(:),y_smoothed(:),'b');
% plot(x(:),y_der_naive(:),'r.');
% plot(x(:),y_derivative(:),'r'); hold off
%
%Test case 2:
% y = [zeros(1,100) 1:100] + 0.5*randn(1,200);
% [y_der,y_smooth] = smooth_derivative(y,10,'linear');
% y_der_naive = [y(:,2)-y(:,1), 0.5*[y(:,3:end)-y(:,1:end-2)], y(:,end)-y(:,end-1)];
% figure(1); plot(y); hold on;  
% plot(y_smooth,'r'); hold off; legend('Data','Smoothed');
% figure(2); plot(y_der_naive(:),'.'); hold on
% plot(y_der,'r');  hold off; legend('Differences','Smooth derivative');

% Wolf 2005

if ~exist('n_window','var'), n_window = 5; end
if ~exist('x','var'), x = 1:size(y,2); end
if ~exist('method','var'), method = 'cubic'; end

y_derivative = nan*ones(size(y));
y_smoothed   = nan*ones(size(y));

for it = 1:size(y,2),
   indices = max(1,it - n_window):min(it+n_window,size(y,2));
   pos = find(it==indices);
   xval = x(indices);
   yval = y(:,indices);
   weights = exp( -((1:length(indices))-pos).^2/(0.5*n_window^2)); 
   weights = weights/sum(weights);
switch method,
  case 'linear',
   G = [ xval' ones(length(indices),1)];
   a = (pinv(diag(weights)*G)*diag(weights)*yval')'; 
   y_derivative(:,it) = a * [ 1         0]';
   y_smoothed(:,it)   = a * [ xval(pos) 1]';
  case 'square',
   G = [ xval'.^2 xval' ones(length(indices),1)];
   a = (pinv(diag(weights)*G)*diag(weights)*yval')'; 
   y_derivative(:,it) = a * [2*xval(pos),    1         0]';
   y_smoothed(:,it)   = a * [  xval(pos).^2, xval(pos) 1]';
  case 'cubic',
   G = [ xval'.^3 xval'.^2 xval' ones(length(indices),1)];
   a = (pinv(diag(weights)*G)*diag(weights)*yval')'; 
   y_derivative(:,it) = a * [3*xval(pos).^2 2*xval(pos),    1         0]';
   y_smoothed(:,it)   = a * [  xval(pos).^3  xval(pos).^2 xval(pos) 1]';
  end
end
