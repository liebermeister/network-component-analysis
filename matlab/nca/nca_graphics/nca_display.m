function nca_display(X,A,B,Xpred,W,netsigns,fontsize,tsf_list,true_sign_fraction)

display('Displaying normalised values');

% mandatory arguments: X, A, B, W

eval(default('netsigns','1','tsf_list','[]','true_sign_fraction','[]','fontsize','3'));
print = 1;

Anorm = A * diag(std(B'));
Bnorm = diag(1./(std(B')+10^-10))*B;

%Anorm = A * diag(max(abs(B')));
%Bnorm = diag(1./max(abs(B')+10^-10))*B;

if netsigns,
  figure(1); clf; set(gca,'FontSize',fontsize);
  im((netsigns+0.2).*W.signs,[],W.gene_names,W.TF_names); colormap('rb_colors')
%  my_xticklabel(1:length(W.TF_names),0.4,W.TF_names,2*fontsize)
  if print, title('Connectivity matrix');   else, noticks; end
  
  %%figure(101);
  %%im(sign(Anorm).*netsigns,[],W.gene_names,W.TF_names);colormap('rb_colors')

  %
  %%if print, set(gca,'FontSize',8); title('Correctness of signs');
  %%else, noticks; end
end

figure(2); clf; set(gca,'FontSize',fontsize);
maxval = 3*nanmedian(nanmedian(abs(Anorm(find(Anorm))))');
im(Anorm,[-maxval,maxval],W.gene_names,W.TF_names); colormap('rb_colors')
%my_xticklabel(1:length(W.TF_names),0.4,W.TF_names,2*fontsize)
if print, set(gca,'FontSize',fontsize); title('Input weight matrix'); %xlabel('TF','FontSize',fontsize);
else, noticks; end

% omit TF that do not change at all:
ind = find(std(Bnorm')>10^-5);
ind2 = sort_by_clustering(Bnorm(ind,:)');
close(10000)
figure(3); clf;
im(Bnorm(ind(ind2),:),[],W.TF_names(ind(ind2)));colormap('rb_colors')

if print, set(gca,'FontSize',fontsize); xlabel('TF activity profiles','FontSize',fontsize);
else, noticks; end

figure(4); clf; 
subplot(1,2,1); im(X,[],W.gene_names); colormap('rb_colors')
if print, set(gca,'FontSize',fontsize); title('Expression data'); else, noticks; end
subplot(1,2,2);
Xpred(~isfinite(X))=nan; 
im(Xpred,[],W.gene_names);colormap('rb_colors')
if print, set(gca,'FontSize',fontsize); title('Expression (reconstructed)'); else, noticks; end

%figure(5)
%hist(tsf_list,15);
%xlabel(sprintf('Percentage of faithfully predicted signs: %s',...
%num2str(true_sign_fraction,2)),'Fontsize',14);
%ylabel('Count numbers under null hypothesis','Fontsize',14);

%figure(6)
%multiple_plot(B,W.TF_names)
