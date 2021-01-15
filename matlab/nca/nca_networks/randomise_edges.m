function WW = randomise_edges(W)

% WW = randomise_edges(W)
% produce a randomised version of a network, preserving degrees of TF and genes

M = W.signs;


%M = double(rand(200,100)<0.1);
%M = sparse(M);

[nr,nc] = size(M); 

for it = 1:10*sum(sum(full(M))),
c = 1+floor(nc*rand(1,2));
y1 = find(M(:,c(1)));
y2 = find(M(:,c(2)));

r1 = y1(1+floor(length(y1)*rand));
r2 = y2(1+floor(length(y2)*rand));

if  sum(sum(M([r1,r2],[c(1),c(2)])-[1 0; 0 1]))==0,
    M([r1,r2],[c(1),c(2)]) = [0 1; 1 0];
end

end

WW.TF_names = W.TF_names;
WW.operon_names = W.operon_names;
WW.signs = full(M);
WW.data = double(WW.signs~=0);
