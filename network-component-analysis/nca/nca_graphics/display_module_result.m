function display_module_result(TFlist,W,X,B,X_pred,Xplus,Xminus,fig1,fig2,names,ffname,linepositions)

printflag = 1;
fitflag = 0;

if ~exist('names','var'), names = 'genes'; end
if ~exist('linepositions','var'), linepositions = []; end

for it =1:length(TFlist),
 
 [Wsub,cg,ct] = pick_subW_TF(W,TFlist{it},0,1);
 figure(fig1+it); 
 clf;
 if iscell(B), 
   thisB = []; for it2=1:length(B), thisB = [thisB {B{it2}(ct,:)}]; end   
   display_curves_and_connectivity(Wsub,X(cg,:),thisB,{...
       Xminus(cg,:),Xplus(cg,:),X_pred(cg,:)},struct('names',names,'linepositions',linepositions,'textsize',12,'textyoffset',0,'textTFxoffset',-2));
 %  im([thisB{1} thisB{2} thisB{3} thisB{4} ]); pause
 else,
   display_curves_and_connectivity(Wsub,X(cg,:),B(ct,:),{...
       Xminus(cg,:),Xplus(cg,:),X_pred(cg,:)},struct('names',names,'linepositions',linepositions,'textsize',12,'textyoffset',0,'textTFxoffset',-2));
 end
 axis off
 name = TFlist{it}{1}; for it2 = 2:length(TFlist{it}), name = [name '/' TFlist{it}{it2}]; end
 fname = TFlist{it}{1}; for it2 = 2:length(TFlist{it}), fname = [fname '-' TFlist{it}{it2}]; end
% title(['Targets of ' name],'FontSize',16);

 if fitflag,
 
 figure(fig2+it),
 clf;
 [ni,nk] = subplot_n(length(cg));
 
 for zz = 1:length(cg),
  subplot(ni,nk,zz);
  it2 = cg(zz);
  relevant = find(isfinite(X(it2,:)));
  x = X_pred(it2,relevant);
  y = X(it2,relevant);
  % FOR UNTRANSFORMED DATA:
%  y = exp(y + log(100))-100;
  ymax_start = 1.5*max(y); km_start = max(x); hc_start = 1;
  kpar = [ymax_start km_start hc_start];
  kparopt = fminsearch(@hill_msr,kpar,optimset('MaxFunEvals',1000),x,y);
  ymax(it) = kparopt(1); km(it) = kparopt(2); hc(it) = kparopt(3);
  y_fit = hill(x,ymax(it),km(it),hc(it));
  X_fit(it,:)= y_fit;
  plot(x,y,'.','color',[0.7,0.7,0.85]); hold on 
  plot(sort(x),sort(y_fit),'r.'); hold off
  noticks
  axis tight
%  title(W.operon_names{it2});
  drawnow
  end
end
%for zz = 1:length(cg),
%  subplot(ni,nk,zz);
%  it2 = cg(zz);
%  relevant = find(isfinite(X(it2,:)));
%  plot(X_pred(it2,relevant),X(it2,relevant),'.','color',[0.7,0.7,0.85]);
%  [d2,order] = sort(X_pred(it2,relevant)+10^-8*rand(1,length(relevant)));
%  d1 = d2(1) + (d2(end)-d2(1))* (0:0.001:1);
%  d3 = interp1(d2, X(it2,relevant(order)),d1,'cubic');
%  hold on; 
%  d4 = conv2(d3,1/300*ones(1,300),'same');
%  plot(d1(150:end-150),d4(150:end-150),'r','linewidth',3);
%  hold off; axis equal;  noticks;  title(W.gene_names{it2});
% end

 if printflag,
  disp( [' Writing image file ' ffname '_NCA_results_TFmodule_' fname '_curves.eps'] )
  print([ffname '_NCA_results_TFmodule_' fname '_curves.eps'],['-f' num2str(fig1+it)],'-depsc');
  if fitflag,
    print([ffname '_NCA_results_TFmodule_' fname '_quality.eps'],['-f' num2str(fig2+it)],'-depsc');
  end
 end
end
