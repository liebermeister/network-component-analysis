function [est,est1,est2] = gfp_data_estimate_error(data,use_time_points,mean_cv);

% [est,est1,est2] = gfp_data_estimate_error(data,use_time_points,mean_cv);


est1 = 0*(10.^(1.5+(0.9/1.5) * (log10(data+1)-2)));  % from 3 gene repeats

est2 = data * diag(mean_cv);
est2=est1;
est3 = 10.^(1 + 0.5 * (log10(data+1)-1));  % from repeated plates
est3(data<10) = 10;
%est3 = repmat(median_dist(use_time_points),size(data,1),1);

est4 = abs(shift_matrix(data,-2)-shift_matrix(data,2));

est = reshape(max([reshape(est4',prod(size(est4)),1)';reshape(est1',prod(size(est1)),1)';reshape(est2',prod(size(est1)),1)' ;reshape(est3',prod(size(est1)),1)']),size(est1,2),size(est1,1))';

%est(est<10)=10;
