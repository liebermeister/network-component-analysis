function display_input_weights(info,W,X,filename,TF_choice,operon_choice,general,lambda_dep,regularise_flag)

% -----------------------------------------------------------
% Show reproducibility of several runs

A = info.A.best;
B = info.B.best;

a           = 1;  % 4 fuer ohne allg. einfluesse, 1 mit
B_undetermined = (W.data(:,a:end)'*isfinite(X)==0);
B           = B .* (1-B_undetermined);

if ~exist('offset','var'), offset=[]; end

if ~isempty(offset),
clear offset
for it=1:2,
  offset(:,D.experiments{it}) = ...
      repmat(nanmean(B(:,D.experiments{it})')',1,length(D.experiments{it}));
end
offset(1:2,:) = 0;
B = B - offset;
end

B = diag(1./sqrt(mean(B.^2,2)))*B;
A = optimise_A_given_B(X,W.data,B,W.signs,0,A,general,lambda_dep,regularise_flag);

[dummy,order] = sort(info.errors);
order = order(1:min(20,length(order)));

% normalise all runs to same mean square B (and signs, if change is possible)
clear Arun Brun Xrun;

for it = 1:length(order); 
  Brun{it} =  info.B.list{order(it)} .* (1-B_undetermined);
  if ~isempty(offset),
    clear offset
    for it2=1:2,
      offset(:,D.experiments{it2}) = ...
	  repmat(nanmean(Brun{it}(:,D.experiments{it2})')',1,length(D.experiments{it2}));
    end
    offset(1:2,:) = 0;
    Brun{it}=Brun{it}-offset; 
  end
  Brun{it} = diag(1./(sqrt(mean(Brun{it}.^2,2))+10^-10))*Brun{it};
  Arun{it} = optimise_A_given_B(X,W.data,Brun{it},W.signs,0,info.A.list{order(it)},general,lambda_dep,regularise_flag);
  [Arun{it},Brun{it}] = nca_adjust_signs(Arun{it},Brun{it},sign(A),'match B',B);
  Xrun{it} = Arun{it} * Brun{it}; Xrun{it} =  Xrun{it} .* isfinite(X);

  Brun_adj{it} = Brun{it};
switch filename, 
  case 'diauxic_casamino__regulonDB',
% for diauxic+casamino
    Brun_adj{it}(4:end,:) =   Brun{it}(4:end,:) -  ( (Brun{it}(4:end,:) - B(4:end,:))  *pinv(Brun{it}(1:3,:)) ) * Brun{it}(1:3,:);
  otherwise,
% for single experiments
    Brun_adj{it}(2:end,:) =   Brun{it}(2:end,:) -  ( (Brun{it}(2:end,:) - B(2:end,:))  *pinv(Brun{it}(1,:)) ) * Brun{it}(1,:);
end
Brun_adj{it} = diag(1./sqrt(mean(Brun_adj{it}.^2,2)))*Brun_adj{it};

end

% ----------------------------------------------------------

clf;
for it=1:length(TF_choice),

  
 TF_index = TF_choice(it);
 
 clear dum dumA dum_adj;
 for it2 = 1:length(order);  
   dum(it2,:)     =  Brun{it2}(TF_index,:);
   dumA(:,it2)    =  Arun{it2}(:,TF_index);   
   dum_adj(it2,:) =  Brun_adj{it2}(TF_index,:);
 end
 d1 = min(dum); d2 = max(dum);
 d1_adj = min(dum_adj); d2_adj = max(dum_adj);
 
subplot(4,1,it);
set(gca,'Fontsize',16);
ind = operon_choice;
relevant = find(W.data(operon_choice,TF_index)~=0);
mediana = median(dumA(ind,:)',1);
m1a  = max(dumA(ind,:)');   m2a  = min(dumA(ind,:)');
m3a  = nanquantile(dumA(ind,:)',0.75);   m4a  = nanquantile(dumA(ind,:)',0.25);
m0a  = A(ind,TF_index)';
if length(relevant),
  h1=  errorbar(relevant, mediana(relevant),mediana(relevant)-m2a(relevant),m1a(relevant)-mediana(relevant),'c.'); hold on;
  h3=  plot(relevant,m0a(relevant),'.','marker','square','color','r','markerfacecolor','r');
  h2=  errorbar(relevant, mediana(relevant),mediana(relevant)-m4a(relevant),m3a(relevant)-mediana(relevant),'b.');hold off;
end
ylabel( W.TF_names{TF_index});
axis([0.5, length(ind)+.5, 1.2*min(m0a) - 0.2*max(m0a), 10^-5+(1.2*max(m0a) -0.2*min(m0a))]);
line( [0.5, length(ind)+.5],[0 0] )
%   noticks
set(gca,'Fontsize',16,'XTick',[]);
end
set(gca,'XTick',1:length(ind),'XTicklabel',W.operon_names(operon_choice));

