function [GFP_p_OD,GFPder_p_OD,OD,GFP] = D_display(D,Wblock,graphics_flag)

% [GFP_p_OD,GFPder_p_OD,OD,GFP] = D_display(D,Wblock,graphics_flag)
% 
% choose 

if ~exist('graphics_flag','var'), graphics_flag = 1; end

if ~exist('Wblock','var'), 
  operon_indices_choice = 1:size(D.GFPder_p_OD,1); 
else,
  operon_indices_choice = find_operon(Wblock.operon_names,D);
end

% ------------------------------
% graphics for the chosen GFP_p_OD

if graphics_flag,
 GFP_p_OD    = normleft_max(D.GFP_p_OD(operon_indices_choice,:));
 GFPder_p_OD = normleft_max(D.GFPder_p_OD(operon_indices_choice,:));
 OD          = normleft_max(D.OD_shifted(operon_indices_choice,:));
 GFP         = normleft_max(D.GFP_shifted(operon_indices_choice,:));

im([GFPder_p_OD GFP_p_OD OD GFP],[],D.gene_names(operon_indices_choice));
 set(gca,'FontSize',8)
% figure(2); im(D.GFP_p_OD,[],Wblock.operon_names,Wblock.TF_names);
% set(gca,'FontSize',8)
 %figure(3); im([D.GFP_p_OD,GFP_p_OD(:,1:50)]);
 %set(gca,'YTick',1:length(Wblock.operon_names),'YTickLabel',Wblock.operon_names,'FontSize',8)
end
