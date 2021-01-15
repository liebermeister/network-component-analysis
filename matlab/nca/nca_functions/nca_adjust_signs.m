function [Anorm,Bnorm] = nca_adjust_signs(A,B,Wblock_signs,method,B_standard)

%function [Anorm,Bnorm] = nca_adjust_signs(A,B,Wblock_signs,method,B_standard)
% (arbitrary yet reproducible) choice of signs for transcription factors
% method 'match B';

if ~exist('method','var'), method = 'match B'; end

switch method,
  case 'from network',
    % signs (by comparison with signs in the network)
    netsigns = zeros(size(Wblock_signs)); 
    netsigns(find(Wblock_signs==1))=1; 
    netsigns(find(Wblock_signs==2))=-1;
    signs = 2*(sign(sum(sign(Anorm).*netsigns))>0)-1;
    scale = ones(size(A,2),1);
    
  case 'positive TF',
    % signs: TF activities should be positive
    signs = sign(sum(Bnorm'));  signs(signs==0)=1;
    scale = ones(size(A,2),1);

  case 'scale to max B',
    signs = ones(size(A,2),1);
    scale = 1./max(abs(B'));

  case 'like liao',    
    % as liao
    signs = sign(diag(B*B_standard')) .* sqrt(mean(B'.^2))';  signs(signs==0)=1;
    val   = sum(A~=0)./sum(abs(A));
    signs = signs.*val;
    scale = ones(size(A,2),1);
    
  otherwise, 
    % TF time courses should match standard - if allowed by network
    scale = 1./(std(B')./std(B_standard'));

    signs = sign(diag(B*B_standard'));
    signs(signs==0)=1;
    signs(find((signs==-1)' .* (sum(abs(Wblock_signs))))) = 1;
%    signs = signs .* sqrt(mean(B'.^2)./mean(B_standard'.^2))';  

   %for it = 1:size(B,1),
   %  factor(it) = median_match(B_standard(it,:),B(it,:)+10^-10);
   %end
   %signs = signs  .* factor';

end
Anorm = A*diag(signs)*diag(1./scale);
Bnorm = diag(scale)*diag(signs)*B;
