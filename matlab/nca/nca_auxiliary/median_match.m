function a = median_match(x,y)

a = exp(fminbnd(@dist,log(min(x./y)),log(max(x./y)),[],x,y));
%a = fminbnd(@dist,nanquantile(x'./y',0.3),nanquantile(x'./y',0.7),[],x,y);

function d = dist(a,x,y)
 d = abs(sum(sign(x-exp(a)*y)));
