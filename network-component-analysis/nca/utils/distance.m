% calculate distance or similarity matrix between columns of v1 and v2 
%
% function d = distance(name,v1,v2)
% name in {'euclidean','correlation','half_correlation',
%          'half_half_correlation','covariance','rank_correlation'}

function d = distance(name,v1,v2)

if ~exist('v2'), v2 = v1; end

s1 = size(v1,2);
s2 = size(v2,2);

switch name

 case 'euclidean'
  if s1==1 & s2==1, d = sqrt(sum((v1-v2).^2));
  else, d = real(sqrt( repmat( sum(v1.^2)',1,s2) - 2 * v1'*v2 + repmat(sum(v2.^2),s1,1) )); end

 case 'evidence_of_coregulation'

  eta=0.1;               % noise level
  nsamples = 200;
  ndim = size(v1,1);
  euc = norm(v1-v2);
  vv  = norm(v2);
  q = randn(ndim,nsamples);
  q = q ./ repmat( sqrt(sum(q.^2)), ndim,1) * vv;
  q = q - repmat(v1,1,nsamples);
  noise = sqrt(2) * eta * randn(1,nsamples);
  q = sqrt(sum(q.^2)) + noise;
  d = sum(q>euc)/nsamples;

  case 'noisy_corrcoef',
  d = noisy_corrcoef([v1 v2],1,0);

 case 'correlation'
  d = corrcoef([v1,v2]);
    if s2 == 1, d = d(1,2); else, d = d(s1+1:end,1:s2); end

 case 'half_correlation'
    if s2 == 1,
    d = corrcoef([v1,v2]);
    d = d(s1+1:end,1:s2);
     d = d * sqrt(sqrt(sum(v2.^2)));
     d = d(1,2);
    else
     d = corrcoef([v1,v2]);
     d = d(s1+1:end,1:s2) * diag(sqrt(sqrt(sum(v2.^2))));
    end

  case 'half_half_correlation'
    if s2 == 1,
        d = corrcoef([v1,v2]);
    d = d(s1+1:end,1:s2);
     d = sqrt(sqrt(sum(v1.^2))) * d * sqrt(sqrt(sum(v2.^2)));
     d = d(1,2);
    else
     d = corrcoef([v1,v2]);
     d = diag(sqrt(sqrt(sum(v2.^2)))) * d(s1+1:end,1:s2) * diag(sqrt(sqrt(sum(v2.^2))));
    end

 case 'covariance'
   d = cov([v1,v2]);
    d = d(s1+1:end,1:s2);

    if s2 == 1, d = d(1,2); end

 case 'rank_correlation'
    [dum,w1] = sort(v1,1);
    [dum,w2] = sort(v2,1);
    d = corrcoef([w1,w2]);
    if s2 == 1, d = d(1,2); else ,     d = d(s1+1:end,1:s2); end

 otherwise 
    d = '';

end
