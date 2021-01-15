function [X,X_error,Xmean,Xplus,Xminus] = make_X_matrices(filename);

% make_X_matrices(network_name,data_name);
% different kinds of preprocessing and network models 
% for amino acid data

file_in1  = ['matched_' filename];
file_in2  = ['identifiable_' filename];
file_out  = ['final_' filename];

cd_data(filename,'analysis');
load(file_in1);
load(file_in2);

display('Computing transformed data matrix X');

%figure(1); W_display(W,'M');
%figure(2); D_display(D);
%figure(3); display_curves_and_connectivity(W,X)

% take log (robust), subtract gen mean

if isfield(D,'GFPder_p_OD_Std_err'),
  data_type = 'transcription_rate';
else
  data_type = 'log_transcript_level';
end

switch data_type,
  
  case 'log_transcript_level',
    X = D.data;
    X_error = ones(size(D.data));
    Xplus   = X+X_error;
    Xminus  = X-X_error;
    Xmean   = X;

  case  'transcription_rate',
    
    if ~isfield(D,'GFPder_p_OD_Std_err'), 
      D.GFPder_p_OD_Std_err = 100 *ones(size(D.GFPder_p_OD));
      disp('Warning: no standard error found. Setting constant value 100');
    end
    
    X = D.GFPder_p_OD;
    X_err = D.GFPder_p_OD_Std_err;

    % fill missing values
    display([' Filling missing values']);
 
    blocks = D.experiments;
    X      = [];
    X_err  = [];
    for it = 1:length(blocks),
      X= [X fill_nan(D.GFPder_p_OD(:,blocks{it}))];
      this_Xerr = D.GFPder_p_OD_Std_err(:,blocks{it});
      this_Xerr(isnan(D.GFPder_p_OD(:,blocks{it}))) = 300;
      X_err = [X_err  this_Xerr];
    end
    
Xplus   = X + X_err;
Xminus  = max(0,X - X_err);
Xminus(isnan(X))=nan;

% log-transformation

X       = log(X+100) - log(100);
Xplus   = log(Xplus+100)  - log(100);
Xminus  = log(Xminus+100) - log(100);
X_error = Xplus-X;
Xmean   = nanmean(X')';

end

cd_data(filename,'analysis');
save(file_out,'X','X_error','Xmean','Xplus','Xminus'); 

display([' Writing file "' file_out '"']);
