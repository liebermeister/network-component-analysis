function Wblock = W_add_general(Wblock,name,sgns)

% Wblock = W_add_general(Wblock,name,signs)
% 
% add a general TF to the network

if ~exist('name','var'), name = 'General input'; end
if ~exist('sgns','var'), sgns = 1;  end

Wblock.TF_names =[{name}; Wblock.TF_names];
Wblock.data = [ones(size(Wblock.data,1),1) Wblock.data];
Wblock.signs = [sgns*ones(size(Wblock.signs,1),1) Wblock.signs];
