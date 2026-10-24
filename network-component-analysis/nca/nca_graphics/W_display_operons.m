function W_display_operons(W,operons,augment)

% W_display_operons(W,operons)
%
% display a subnetwork of the network described in W containing 
% only the operons given (if augment = 0, default), or the operons given
% and all operons targeted by their transcription factors (if augment ~=0)

if ~exist('augment','var'), augment = 0; end

switch augment,
  case    0, Wsub = pick_subW_operons(W,operons,0);
  otherwise, Wsub = pick_subW_TF(W,W.TF_names(find(sum(W.data(find_operon(operons,W),:)))),0,0);
end
   
W_display(Wsub,'network',struct('xTF',-0.1,'xOP',1.01);
