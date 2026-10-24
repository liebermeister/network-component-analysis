function matrix=num2bit(id,n)

if nargin<2,
  n=max(1,ceil(log2(max(id))));
end

for it = 1:length(id),

b = [];
for i=0:n-1,
   b=[b bitand(id(it),2^(i))];
end;
b = b>0;

matrix(it,:)=b;

end
