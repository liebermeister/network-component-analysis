function p = traveling_salesman(data,nit,method,seed)

% p = traveling_salesman(data,nit,method,seed)
%
% traveling salesman algorithm from /usr/global/matlab6/toolbox/matlab/demos/travel.m
% points: COLUMNS of data
%
% test
% x=rand(2,30); p=traveling_salesman(x); plot(x(1,p),x(2,p))

% -----------------------------------------
% my own changes:

eval(default('nit','1000','seed','0','method','''euclidean'''));

rand('state',seed)

switch method,
  case 'euclidean',
    distmatrix = sqrt(repmat(sum(data.^2),size(data,2),1) + repmat(sum(data.^2),size(data,2),1)' - data'*data);
  case 'distance',
    distmatrix = data;
end

% -----------------------------------------
% original code:

npts = size(data,2);   

% Generate an initial random path between those cities

%  p = randperm(npts);
p = 1:npts;

 for it=1:nit
   len=LocalPathLength(p,distmatrix);
   lenhist=len;
   
   T = 10*len*nit/it;
     
% Try a point for point swap
% ========================
   swpt1=floor(npts*rand)+1;
   swpt2=floor(npts*rand)+1;
   
   swptlo=min(swpt1,swpt2);
   swpthi=max(swpt1,swpt2);
   
   order=1:npts;
   order(swptlo:swpthi)=order(swpthi:-1:swptlo);
      pnew = p(order);
      
      lennew=LocalPathLength(pnew,distmatrix);
      if lennew<len | rand<exp(-T*(lennew-len)),
         p=pnew;
         len=lennew;
         drawFlag=1;
      end;
      % ========================
      
      % Try a single point insertion
      % ========================
      swpt1=floor(npts*rand)+1;
      swpt2=floor((npts-1)*rand)+1;
      
      order=1:npts;
      order(swpt1)=[];
      order=[order(1:swpt2) swpt1 order((swpt2+1):(npts-1))];
      pnew = p(order);
      
      lennew=LocalPathLength(pnew,distmatrix);
      if lennew<len | rand<exp(-T*(lennew-len)),
         p=pnew;
         len=lennew;
      end

      
      % Try a loop inversion
      % ========================
      swpt1=floor(npts*rand)+1;
      swpt2=swpt1+floor((npts-swpt1)*rand);
      pnew = p;
      pnew(swpt1:swpt2)=fliplr(p(swpt1:swpt2));
      lennew=LocalPathLength(pnew,distmatrix);
      if lennew<len | rand<exp(-T*(lennew-len)),
         p=pnew;
         len=lennew;
      end

end

% ----------------------------------------------------

function total=LocalPathLength(p,distmatrix);

% Calculate current path length for traveling salesman problem.
% This function calculates the total length of the current path
% p in the traveling salesman problem.

npts = size(p,2);

% This is a vectorized distance calculation
%
% We're creating two vectors: p and p([end 1:(end-1)]
% The first is the list of first cities in any leg, and the second
% is the list of destination cities for that same leg.
% (the second vector is an element-shift-right version of the first)
% 
% We then use column indexing into distmatrix to create a vector of
% lengths of each leg which can then be summed.

total=sum(distmatrix([(p-1)*npts + p([end 1:(end-1)])]));

