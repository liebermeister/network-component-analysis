%DENSITYPLOT scatterplot using gaussian kernel density estimator
%
%   function densityPlot(X,limits,bins,sigma)
%
%   arguments:
%
%     X       2*n data matrix
%     limits  [xmin xmax ymin ymax]    
%     bins    [nx ny] numbers of pixels
%     sigma   width (standard deviation) of gaussian kernel
%
%   (c) Wolfram Liebermeister 2001

function densityPlot(X,limits,bins,sigma)

xmin=limits(1);
xmax=limits(2);
ymin=limits(3);
ymax=limits(4);

ndata = size(X,2);

nx=bins(1);
ny=bins(2);

bx=(xmax-xmin)/nx;
by=(ymax-ymin)/ny;

nbx = ceil(3*sigma/bx);
nby = ceil(3*sigma/by);

kx = bx*(-nbx:nbx);
ky = by*(-nby:nby);

kernel=exp(-(kx/sigma).^2)'*exp(-(ky/sigma).^2);
kernel=kernel./sum(sum(kernel));

%figure(1)
%plot( X(1,:) , X(2,:),'.');
%axis([xmin xmax ymin ymax]);

XX  = ceil( [nx/(xmax-xmin) 0; 0 nx/(ymax-ymin)] * ...
            (X-repmat([xmin; ymin],1,ndata)));

XX=XX(: , find( (min(XX)>0) .* (XX(1,:)<=nx) .* (XX(2,:)<=ny) ) );

matrix=zeros(nx,ny);

for k=1:size(XX,2);
%XX(1,k)
%     XX(2,k)
  matrix( XX(1,k),XX(2,k) ) = matrix( XX(1,k),XX(2,k) )+1;
end

matrix=flipud(matrix');

matrix=conv2(matrix,kernel','same');

imagesc(log(0.1+matrix));
set(gca,'XTick',[1 nx/2 nx]);
set(gca,'XTicklabel',[xmin xmin+(xmax-xmin)/2 xmax]);
set(gca,'YTick',[1 ny/2 ny]);
set(gca,'YTicklabel',[ymin ymin+(ymax-ymin)/2 ymax]);
