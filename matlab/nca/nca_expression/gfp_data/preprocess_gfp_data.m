%function [DD,info] = preprocess_gfp_data(D,parameters)
%
%preprocess time curves from gfp reporter experiments
%
%INPUTS
% D  structure with fields
% .gene_names         column list of strings
% .GFP                matrix of gfp fluorescence values, each gene is a row
% .OD                 matrix of optical density values, each gene is a row
% .gene_annotations   (optional field) column list of strings 
%
% parameters (list of parameter names, each followed by the respective value):
%  'graphics_flag'       flag: show graphics? (default 0 =no)
%  'remove_genes_flag'   flag: removes genes with suspicious data? (default 1 =yes)
%  'use_time_points'     vector, indices of time points that should be used 
%                        (default: all timepoints starting from point 21)
%  'max_shift'           maximal time shift allowed for OD alignment (default 20)
%  'smoothing_window'    window size for cubic smoothing of time curves (default 20)
%  'unique_genes'        flag: should genes with the same name be averaged over? (default 1 =yes)
%  'tmin', 'tmax'        indices of data points between which the curves are to be matched
%                        (indices within the original data set, not use_time_points, default 10,50)
%  'match_what'          'growth_rate','OD' (default 'growth_rate');         
%
%OUTPUTS
% DD  structure containing the preprocessed data with fields
% .gene_names          gene names as provided by the input
% .GFP_p_OD            shifted + smoothed GFP/OD
% .GFPder_p_OD         shifted + smoothed (dGFP/dt)/OD
% .GFPder_p_OD_Std_err (rough!!) estimate of the errors in .GFPder_p_OD
% .GFP_shifted         shifted version of GFP curves
% .OD_shifted          shifted version of OD curves
% .growth_rate         shifted version of (dOD/dt)/OD curves
% .gene_annotations    (if no annotations were provided -> contains the gene names)
% .original_indices    vector of gene indices in the original data
% .shifts              vector of time shifts applied to the curves
%
% info  structure containing additional information with fields
%  .no_increasing_OD     index vector denoting genes with decreasing OD
%  .no_small_shift           index vector denoting genes with extreme required time shift
%  .no_increasing_GFP    index vector denoting genes with decreasing GFP
%  .no_considerable_gfp  index vector denoting genes with low GFP
%  .mean_OD              mean OD curve
%  .scale                scaling factors for OD curves (for the alignment)
%
%USAGE EXAMPLE
% D.GFP        = gfp;
% D.OD         = od;
% D.gene_names = gene_names;
% [DD,info] = preprocess_gfp_data(D);
%% or, with some parameters,
% [DD,info] = preprocess_gfp_data(D,'graphics_flag',0,'use_time_points',21:size(D.OD,2));
%
%WHAT THE SCRIPT ACTUALLY DOES:
% 1. determine a time shift for each OD curve, such that, in combination with an optimal
%    rescaling (determined seperately for each curve), the curve matches optimally the 
%    median of all curves between tmin (default 10) and tmax (default 50)
% 2. Determine suspicious genes fulfilling at least one out of the four criteria mentioned
%    in the explanation of .info above
% 3. Remove them if this is indicated by the flag 'remove_genes'
% 4. Compute smoothed versions and smoothed derivatives of OD and GFP curves (by using 
%    'smooth_derivative.m'. Shift them according to the time shift determined in (1.)
%    and compute GFP/OD and (dGFP/t)/OD
% 5. Compute error estimates of (dGFP/t)/OD (see the script 'gfp_data_estimate_error')
% 6. Average over repeated occurrences of genes (see the script 'gfp_data_unify_genes')
%    if this is indicated by the flag 'unique_genes'

% wolf, 2005

function [DD,info] = preprocess_gfp_data(D, s1, v1, s2, v2, s3, v3, ...
      s4, v4, s5, v5, s6, v6, s7, v7, s8, v8, ...
      s9, v9, s10, v10, s11, v11, s12, v12, ...
      s13, v13, s14, v14, s15, v15)


% ------------------------------------
% set parameters

graphics_flag     = 0;
remove_genes_flag = 1;
use_time_points   = 21:size(D.GFP,2);
max_shift         = 20;
smoothing_window  = 20;
unique_genes      = 1;
n_wells           = 96;
n_plates          = size(D.GFP,1)/n_wells;
match_what        = 'growth_rate';
tmin              = 10;
tmax              = 50;
consider_plate_mean_shift = 0;    
match_what        = 'growth_rate';

% actual shift = mean of gene shift and plate shift. currently not used:

if(rem(nargin-1,2)==1)
  error('Optional parameters should always go by pairs');
else,
  for i=1:(nargin-1)/2
    % get the name and value of parameter
    str_param = eval (['s' int2str(i)]);
    val_param = eval (['v' int2str(i)]);
    switch str_param, 
      case 'graphics_flag', graphics_flag = val_param;
      case 'remove_genes_flag', remove_genes_flag = val_param;
      case 'use_time_points', use_time_points = val_param;
      case 'smoothing_window', smoothing_window = val_param;
      case 'max_shift', max_shift = val_param;
      case 'unique_genes', unique_genes = val_param;
      case 'tmin', tmin = val_param;
      case 'tmax', tmax = val_param;
      case 'match_what', match_what = val_param;
    end
  end
end

% select the time window indicated in use_time_points

OD  = D.OD(:,use_time_points);
GFP = D.GFP(:,use_time_points);
tmin = sum(use_time_points<=tmin);
tmax = sum(use_time_points<=tmax);

if ~isfield(D,'gene_names'), D.gene_names=numbered_names('gene_',size(GFP,1)); end

% ------------------------------------------------------
% determine the curves to be matched

switch match_what,
  
  case 'growth_rate',
    [ODder,OD_smoothed]   = smooth_derivative(OD,smoothing_window,'cubic');
    ODder_p_OD_smoothed  = ODder ./ OD;
    growth_rate          = my_bfilt(fill_nan(ODder_p_OD_smoothed));
    growth_rate(growth_rate<0) = 0;
    match_curves = growth_rate;

  case 'OD',
    match_curves = OD;    
end
 
mean_curve   = median(match_curves);


% ------------------------------------------------------
% shift and rescale gfp and od curves

fprintf('Computing the curve shifts\n');

clear OD_shifted GFP_shifted shift scale

for it = 1:size(OD,1);
 [shift(it,1),scale(it)] =  ...
     transform_OD_timecourse(match_curves(it,:)',mean_curve',[],max_shift,tmin,tmax);
 OD_shifted(it,:)  = shift_matrix (OD(it,:),shift(it));
 GFP_shifted(it,:) = shift_matrix(GFP(it,:),shift(it));
end


% -----------------------------------------
% find and remove invalid genes

for it = 1:size(GFP,1),
  [ dum,order] = sort(GFP(it,:)); %  plot(order,'.');
  c(it) = my_corrcoef(1:length(order),order);
end

increasing_OD    = (OD(:,end)'./max(OD')>0.9);
small_shift      = (abs(shift)<max_shift); 
increasing_GFP   = (c > 0.8);
considerable_gfp = (median(GFP') > nanquantile(median(GFP')',0.09));

ind_valid = find(  increasing_OD .*  small_shift' .*  increasing_GFP .*  considerable_gfp);
ind_invalid = setdiff(1:size(OD,1),ind_valid);

fprintf('The following genes look suspicious:\n');
table(sort(D.gene_names(ind_invalid))',0)

fprintf('because...\n... OD decreases:\n');
table(sort(D.gene_names(find(1-increasing_OD)))',0)
fprintf('\n... requested curve shift is too large:\n');
table(sort(D.gene_names(find(1-small_shift)))',0)
fprintf('\n... GFP doesnt increase:\n');
table(sort(D.gene_names(find(1-increasing_GFP)))',0)
fprintf('\n... GFP is too low:\n');
table(sort(D.gene_names(find(1-considerable_gfp)))',0)


% ----------------------------------------------------------
% use mean shifts per plate for a more robust estimation of the shift 
% (currently not used)

if consider_plate_mean_shift,
  shift_uncorr       = shift;
  OD_shifted_uncorr  = OD_shifted;
  GFP_shifted_uncorr =  GFP_shifted;
  shift = shift_uncorr; shift(ind_invalid) = nan;
  dum = reshape(shift,n_wells,n_plates); % figure(111); im(dum);
  dum = round(0.5*dum + 0.5*repmat(nanmean(dum),size(dum,1),1)); % figure(112); im(dum);
  dum = reshape(dum,size(shift));
  shift(ind_valid ) = dum(ind_valid);
  shift(ind_invalid ) = shift_uncorr(ind_invalid);
  clear OD_shifted_corr GFP_shifted_corr
  for it = 1:size(OD,1);
    OD_shifted(it,:) = shift_matrix(OD(it,:),shift(it));
    GFP_shifted(it,:)= shift_matrix(GFP(it,:),shift(it));
  end
end


% -------------------------------------------------------

fprintf('\nComputing smooth derivatives gfp(t)\n');

if remove_genes_flag,
  keep_genes = ind_valid;
else,
  keep_genes = 1:length(D.gene_names);
end 

[GFPder,GFP_smoothed] = smooth_derivative(GFP(keep_genes,:),smoothing_window,'cubic');
[ODder,OD_smoothed]   = smooth_derivative(OD(keep_genes,:),smoothing_window,'cubic');

GFP_p_OD_smoothed = GFP_smoothed ./ OD(keep_genes,:);
GFP_p_OD_smoothed = my_bfilt(fill_nan(GFP_p_OD_smoothed));

GFPder_p_OD       = GFPder ./ OD(keep_genes,:);
GFPder_p_OD       = my_bfilt(fill_nan(GFPder_p_OD));
GFPder_p_OD(find(GFPder_p_OD<=0))=0;

growth_rate       = ODder ./ OD(keep_genes,:);
growth_rate       = my_bfilt(fill_nan( growth_rate));
growth_rate(growth_rate<0) = 0; 
 
for it = 1:length(keep_genes);
  GFP_p_OD_smoothed(it,:) = shift_matrix(GFP_p_OD_smoothed(it,:),shift(it));
  GFPder_p_OD(it,:)       = shift_matrix(GFPder_p_OD(it,:),shift(it));
  growth_rate(it,:)       = shift_matrix(growth_rate(it,:),shift(it));
end


% ----------------------------------------------------------
% assembling the output variables DD and info


clear DD;
DD.GFP_shifted        = GFP_shifted(keep_genes,:);
DD.OD_shifted         = OD_shifted(keep_genes,:);
DD.GFP_p_OD           = GFP_p_OD_smoothed;
DD.GFPder_p_OD        = GFPder_p_OD;
DD.growth_rate        = growth_rate;
DD.gene_names         = D.gene_names(keep_genes);
if isfield(D,'gene_annotations'), DD.gene_annotations   = D.gene_annotations(keep_genes);
else, DD.gene_annotations   = D.gene_names(keep_genes); end
DD.shifts             = shift(keep_genes);

if remove_genes_flag,
 DD.original_indices   = keep_genes';
else,
 DD.original_indices = 1:size(GFP_shifted,1);
end 

mean_cv = gfp_data_determine_CV(DD);

if unique_genes, DD = gfp_data_unify_genes(DD); end

[DD.GFPder_p_OD_Std_err,est1,est2] = gfp_data_estimate_error(DD.GFPder_p_OD,use_time_points,mean_cv);

info.no_increasing_OD    = find(increasing_OD==0);
info.no_small_shift      = find(small_shift==0);
info.no_increasing_GFP   = find(increasing_GFP==0);
info.no_considerable_gfp = find(considerable_gfp==0);
info.mean_OD             = mean_curve;
info.scale               = scale;

% -------------------------
% THE REST IS ONLY GRAPHICS
% -------------------------





% ---------------------------------------------------------
% statistics over plates and positions


if graphics_flag,
figure(1)
plot(OD'); ;%set(gca,'Yscale','Log');
%mean_curve = exp(median(log(OD)));
hold on; plot(mean_curve,'yo'); hold off

figure(2)
plot((diag(scale)*OD_shifted_uncorr)'); hold on; 
plot( mean_curve,'y.');  hold off

figure(4);
subplot(2,2,1); im(1-reshape(increasing_OD,n_wells,n_plates)); 
xlabel('plates'); ylabel('positions'); title('Decreasing OD?');
subplot(2,2,2); im(1-reshape(small_shift,n_wells,n_plates)); 
xlabel('plates'); ylabel('positions'); title('Big shifts?');
subplot(2,2,3); im(1-reshape(increasing_GFP,n_wells,n_plates)); 
xlabel('plates'); ylabel('positions'); title('Decreasing GFP?');
subplot(2,2,4); im(1-reshape(considerable_gfp,n_wells,n_plates)); 
xlabel('plates'); ylabel('positions'); title('GFP above background?');

%figure(105); clf; im(reshape(shift_uncorr,n_wells,n_plates)); xlabel('plates'); ylabel('positions');
%figure(106); clf; im(reshape(log(scale),n_wells,n_plates),[-1,1]); xlabel('plates'); ylabel('positions');
%for repeats: scales & shifts are related, scales are robust agains repeat, shifts less robust
%figure(108); clf; plot(log(scale(:)),shift_uncorr(:),'.');
%figure(109);  clf; plot(log(scale(1:96)),log(scale(97:end)),'.');
%figure(110); clf; plot(shift_uncorr(1:96),shift_uncorr(97:end),'.');

figure(5);
dum = reshape(shift_uncorr,n_wells,n_plates);
clf; for it = 1:n_plates,  subplot(4,6,it); im(reshape(dum(:,it),12,8),[-40,40]); colorbar off; end

figure(6);
dum = reshape(log(scale),n_wells,n_plates);
clf; for it = 1:n_plates,  subplot(4,6,it); im(reshape(dum(:,it),12,8),[-1,1]); colorbar off; end

figure(7)
plot(OD_shifted(ind_valid,:)'); hold on; 
plot(mean_curve,'y.'); 
plot(nanmedian(OD_shifted(ind_invalid,:)),'m.');  hold off

figure(8)
plot((diag(scale(ind_valid))*(fill_nan(GFP_shifted(ind_valid,:))))');
plot((diag(scale(ind_valid))*(fill_nan(OD_shifted(ind_valid,:))))'); set(gca,'YScale','log');

figure(9)
plot((diag(scale(ind_valid))*(fill_nan(GFP_shifted_uncorr(ind_valid,:))))');
plot((diag(scale(ind_valid))*(fill_nan(OD_shifted_uncorr(ind_valid,:))))'); set(gca,'YScale','log');

figure(10);
 %im(GFPder_p_OD(1:100,:)); 
im(GFPder_p_OD); 

end

% plot(median(GFP)); 
% plot(sort(median(GFP'))); set(gca,'YScale','log');
% nanquantile(median(GFP')',0.09)
