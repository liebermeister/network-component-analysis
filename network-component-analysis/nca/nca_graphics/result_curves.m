function result_curves(Wsub_block,operons_sub,TF_sub,X,A,TF_act,Xpred_sub);

% result_curves(Wsub_block,operons_sub,TF_sub,X,A,TF_act,Xpred_sub);

TF_names_sub     = Wsub_block.TF_names;
operon_names_sub = Wsub_block.operon_names;
Xsub      = X(operons_sub,:);
TF_actsub = TF_act(TF_sub,:);
Xredsub   = A(operons_sub,TF_sub)*TF_actsub;

figure(8); close(8)
figure(8)

[pos_TF,pos_operon] = W_display(Wsub_block); hold on;

xvalues = - 1.2 + 1 * (0:1/(size(TF_act,2)-1):1);
ysize = 0.1*(max(pos_operon)-1);
for it = 1:length(Wsub_block.TF_names),
   plot( xvalues, pos_TF(it) + ysize*TF_actsub(it,:),'r');
end

xvalues = 1.3 + 1 * (0:1/(size(TF_act,2)-1):1);
ysize = 0.9;
for it = 1:length(Wsub_block.operon_names),
  plot(xvalues, pos_operon(it) + ysize*Xsub(it,:),'b');
  plot(xvalues(1), pos_operon(it) + ysize*Xsub(it,1),'b.');
  plot(xvalues, pos_operon(it) + ysize*Xredsub(it,:),'r');
if exist('Xpred_sub','var'),
  plot(xvalues, pos_operon(it) + ysize*Xpred_sub(it,:),'m');
end
end

axis tight
hold off;

return

figure(9);
[ni,nk] = subplot_n(size(Xsub,1));
 for it=1:size(Xsub,1),
  set(gca,'Fontsize',16);
  subplot(ni,nk,it);
  plot(Xsub(it,:),'r'); hold on;
  plot(Xsub(it,:),'r.');
  plot(Xpred_sub(it,:),'b');
  plot(X_pred_onlysub(it,:),'--g'); hold off;
  ylabel(operon_names_sub{it});
						noticks;
end
%legend('Original','reconstructed','restricted');
%title('(d/dt GFP)/OD')


figure(10);
[ni,nk]=subplot_n(length(TF_sub));
for it=1:length(TF_sub),
  subplot(ni,nk,it);
  set(gca,'Fontsize',16);
  plot(TF_actsub(it,:));
  ylabel(TF_names_sub{it});
%  title('Reconstructed TF activities')
noticks
end
