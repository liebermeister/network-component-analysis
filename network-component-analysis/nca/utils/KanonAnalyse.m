% function [profil1,profil2,eival1,kanvar1,kanvar2] = KanonAnalyse(X1,X2,nkanvar)
% canonical correlation analysis (using transposed matrices)
% requires matrices X1,X2 with individuals as row vectors

function [profil1,profil2,eival1,kanvar1,kanvar2] = KanonAnalyse(X1,X2,nkanvar)

ndata = size(X1,1);
X1 = X1 - repmat(mean(X1,1),ndata,1);
X2 = X2 - repmat(mean(X2,1),ndata,1);

if ~exist('nkanvar','var'),
  nkanvar = min(size(X1,2),size(X2,2));
end

V11  = X1' * X1;
V12  = X1' * X2;
V21  = V12';
V22  = X2' * X2;
iV11 = pinv(V11);
iV22 = pinv(V22);
Sigma1 = iV11 * V12 * iV22 * V21;
Sigma2 = iV22 * V21 * iV11 * V12;

[eivec1,eival1]=eig(Sigma1);
[eival1,index1]=sort(diag(eival1));
eival1 = flipud(eival1);
eivec1 = eivec1(:,flipud(index1));
eival1 = eival1(1:nkanvar);
eivec1 = eivec1(:,1:nkanvar);

[eivec2,eival2]=eig(Sigma2);
[eival2,index2]=sort(diag(eival2));
eival2 = flipud(eival2);
eivec2 = eivec2(:,flipud(index2));
eival2 = eival2(nkanvar);
eivec2 = eivec2(:,1:nkanvar);

kanvar1 = X1 * eivec1;
kanvar2 = X2 * eivec2;

profil1 = pinv(eivec1)';
profil2 = pinv(eivec2)';

korrVorzeichen=zeros(1,nkanvar);
for k=1:nkanvar
  korrVorzeichen(k)= sign( (kanvar1(:,k)-mean(kanvar1(:,k)))'*(kanvar2(:,k)-mean(kanvar2(:,k))) );
end

profil2 = profil2 * diag( korrVorzeichen );
kanvar2 = kanvar2 * diag( korrVorzeichen );
