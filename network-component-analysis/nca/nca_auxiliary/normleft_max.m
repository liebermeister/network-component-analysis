function M=normleft_max(M);

for it = 1:size(M,1),
  M(it,:) = 1./max(abs(M(it,:))) * M(it,:);
end
