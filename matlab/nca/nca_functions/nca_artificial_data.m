function [W_true,A_true,B_true,X,X_true] = nca_artificial_data(n_TF,n_genes,n_timepoints,ptot,sigma)

% [W_true,A_true,B_true,X,X_true] = nca_artificial_data(n_TF,n_genes,n_timepoints,ptot,sigma)
%
% make artificial data to test nca algorithms
% no arguments are required

if ~exist('n_TF','var'),         n_TF    = 10; end  % use even number (simulation of B_true)
if ~exist('n_genes','var'),      n_genes = 30; end
if ~exist('n_timepoints','var'), n_timepoints = 200; end
if ~exist('ptot','var'),         ptot     = 0.5;end  % percentage of connections with nonzero prior
if ~exist('sigma','var'),        sigma    = 0.2; end

% -----------------
% P_chi: all potential connections
% P probabilities of connections
% W actual connections

stop = 0;
while ~stop;
 P_chi = rand(n_genes,n_TF)<=ptot;
 P = P_chi.*rand(size(P_chi));
 W_true.data = double(rand(size(P))<P);
 if (sum(sum(W_true.data)~=0) == n_TF) & (sum(sum(W_true.data')~=0) == n_genes), stop = 1; end
end


W_true.signs = W_true.data .* sign(randn(size(W_true.data)));

W_true.TF_names   =cellstr( [repmat('TF ',n_TF,1) num2str((1:n_TF)')]);
W_true.gene_names   =cellstr( [repmat('gene ',n_genes,1) num2str((1:n_genes)')]);
W_true.operon_names   =cellstr( [repmat('operon ',n_genes,1) num2str((1:n_genes)')]);
W_true.operon_abbr   =cellstr( [repmat('operon ',n_genes,1) num2str((1:n_genes)')]);
W_true.operon_flags   = ones(n_genes,1);

W_true = make_W_identifiable(W_true);

W_true.data(:,1) = 1;
W_true.signs(:,1) = 1;

[n_genes,n_TF] = size(W_true.data);

% -----------------

A_true = W_true.signs.*abs(randn(size(W_true.data)));

%B_true = [cos( (5*(1:n_TF/2)/(n_TF/2))' *  2*pi*(0:n_timepoints-1)/(n_timepoints-1) );
%          sin( (5*(1:n_TF/2)/(n_TF/2))' *  2*pi*(0:n_timepoints-1)/(n_timepoints-1) )];

mn = (0:(n_timepoints/(n_TF-1)):n_timepoints)';
st = n_timepoints/n_TF;

B_true = exp( - (repmat( (1:n_timepoints),n_TF,1) - repmat(mn,1,n_timepoints)).^2 / (2*st.^2) );

X_true = A_true*B_true;
X = A_true*B_true+sigma*randn(n_genes,n_timepoints);

%figure(1);
%subplot(1,3,1); im(X); subplot(1,3,2); im(A_true); subplot(1,3,3); im(B_true); 
