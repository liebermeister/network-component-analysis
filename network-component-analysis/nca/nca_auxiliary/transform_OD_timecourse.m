function [shift_opt,scale_opt,x_transformed,gfp_transformed] = transform_OD_timecourse(x,x_std,gfp,max_shift,tmin,tmax);

% [shift_opt,scale_opt,x_transformed,gfp_transformed] = transform_OD_timecourse(x,x_std,gfp,max_shift);

z=1;

if ~exist('gfp','var'); gfp = []; end
if ~exist('max_shift','var'), max_shift = 50; end
if ~exist('tmin','var'), tmin = 1; end
if ~exist('tmax','var'), tmax = length(x); end

for shift =-max_shift:max_shift,

sh(z)=shift;

%if shift > 0,
%    a = x(1:end-shift);
%    b = x_std(shift+1:end);
%end

%if shift < 0,
%  a =  x(-shift+1:end);
%  b = x_std(1:end+shift);
%end

if shift > 0,
    a = x(max(tmin-shift,1):tmax-shift);
    b = x_std(tmin+max(shift-tmin+1,0):tmax);
end

if shift < 0,
  a =  x(tmin-shift:min(tmax-shift,length(x)));
  b = x_std(tmin : tmax-(max(tmax-shift-length(x),0)) );
end

if shift == 0,
  a = x(tmin:tmax);
  b = x_std(tmin:tmax);
end

valid = isfinite(a).*isfinite(b);
a = a(find(valid));
b = b(find(valid));

scale(z) = (b'*b)/(a'*b);
 loss(z) = mean((scale(z)*a-b).^2);

z=z+1;
end

[dum,index]=min(loss);

shift_opt = sh(index);
%shift_opt=0;
scale_opt = scale(index);

if nargout>2,

x_transformed = x*scale_opt;

if ~isempty(gfp); gfp_transformed = gfp*scale_opt; end

if shift_opt < 0,
    x_transformed = x_transformed(1-shift_opt:end); x_transformed(length(x)+shift_opt+1:length(x))=nan;
if exist('gfp_transformed','var');    gfp_transformed = gfp_transformed(1-shift_opt:end); gfp_transformed(length(gfp)+shift_opt+1:length(gfp))=nan; end
end

if shift_opt > 0,
  x_transformed(shift_opt+1:end) =  x_transformed(1:end-shift_opt); x_transformed(1:shift_opt)=nan;
  if exist('gfp_transformed','var');  gfp_transformed(shift_opt+1:end) =  gfp_transformed(1:end-shift_opt); gfp_transformed(1:shift_opt)=nan; 
  end
end
end
