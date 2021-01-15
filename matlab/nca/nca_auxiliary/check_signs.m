function [netsigns,total_elements,true_sign_fraction,mem] = check_signs(Wblock,A)

l  = label_names(unique(Wblock.operon_names),Wblock.operon_names);
l2 = setdiff(1:size(Wblock.data,2),find_TF('General input',Wblock));
%l2 = 1:size(Wblock.data,2);

fprintf('Permutation test for prediction of weight signs\n');

netsigns = zeros(size(Wblock.signs)); 
netsigns(find(Wblock.signs==1))=1; 
netsigns(find(Wblock.signs==2))=-1;
netsigns_red = netsigns(l,l2);

total_elements=sum(sum(abs(netsigns_red)));

mem=[];
for it =1:500;

netsigns = zeros(size(Wblock.signs)); 
netsigns(find(Wblock.signs==1))=1; 
netsigns(find(Wblock.signs==2))=-1;
netsigns_red = netsigns(l,l2);

Asigns_red = A(l,l2);

 dum = netsigns_red(find(netsigns_red));
 netsigns_red(find(netsigns_red)) = dum(randperm(length(dum)));
 signs = 2*(sum(Asigns_red.*netsigns_red)>=0)-1;
 Asigns_red = Asigns_red*diag(signs);
 true_sign_fraction =  (total_elements + Asigns_red(:)'*netsigns_red(:))/(2*total_elements);
 mem=[mem true_sign_fraction];
end

%hist(mem,20)

netsigns = zeros(size(Wblock.signs)); 
netsigns(find(Wblock.signs==1))=1; 
netsigns(find(Wblock.signs==2))=-1;
netsigns_red = netsigns(l,l2);
%im(netsigns_red)

Asigns_red = A(l,l2);
%im(Asigns_red)

 dum = netsigns_red(find(netsigns_red));
% netsigns_red(find(netsigns_red)) = dum(randperm(length(dum)));
 signs = 2*(sum(Asigns_red.*netsigns_red)>=0)-1;
 Asigns_red = Asigns_red*diag(signs);
 true_sign_fraction =  (total_elements + Asigns_red(:)'*netsigns_red(:))/(2*total_elements);


