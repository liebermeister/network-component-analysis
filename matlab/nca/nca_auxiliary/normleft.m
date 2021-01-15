function M=normleft(M);
% M=normleft(M);

for it = 1:size(M,1),
  m(it) = 1./(10^-10+sqrt(nanmean(M(it,:).^2,2)));
end

M = M .* repmat((m'),1,size(M,2));
