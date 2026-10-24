function p = local_minimum(x)

diff = x(2:end) - x(1:end-1);

p = find( sign(diff(2:end)) + sign(-diff(1:end-1))==2)+1;