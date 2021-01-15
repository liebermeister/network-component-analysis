function [pos_TF,pos_operon] = W_display(W,method,options)

% [pos_TF,pos_operon] = W_display(W,method,options)
%
% Display a transcriptional network
%
% if method = 'N' -> matrix representation (including signs)
% if method = 'M' -> matrix representation (no signs)
% if method = 'R' -> matrix representation (no signs) without names
% otherwise: network representation
%
% fields of 'options' contain the graphics parameters
%
% 'arrow_numbers',0
% 'textsize',11
% 'blocks',1
% 'xTF',0  %  arrows x value left
% 'xOP',1  %  arrows x value right ; % the y values range from 1 to # operons
% 'net_width', 1 overrides xOP: xOP is set to xTF + net_width
% 'lw',1
% 'arrow_values',[]
% 'names','genes'
% 'textyoffset',0
% 'textTFxoffset',0
% 'textOPxoffset',0;
% 'position' [0.1 0.1 0.85 0.8]

cc   = rb_colors;
red  = cc(1,:);
blue = cc(end,:);
grey = [0.85 0.85 0.9];

show_A = 0;

target_names = W.gene_names;

if ~exist('method','var'), method = ''; end

default = struct('arrow_number',0,'textsize',nan,'blocks',0,'xTF',0,'xOP',1,'lw',1,'arrow_values',[],'names','genes','textyoffset',0,'textTFxoffset',0,'textOPxoffset',0,'position', [0.1 0.1 0.85 0.8],'net_width',1,'A_values',[],'normalise_A',0,'no_gene_labels',0);

if ~exist('options','var'), options = struct; end

ff = fields(default);
for it = 1:length(ff),
  if ~isfield(options,ff{it}),
    options=setfield(options,ff{it},getfield(default,ff{it}));
  end
end

if length(options.net_width),
  options.xOP = options.xTF + options.net_width;
end

if isnan(options.textsize), options.textsize = 1/length(W.gene_names); end

switch options.names,
  case 'operons',
    target_names = W.operon_names;
  case 'both'
   for it=1:length(W.operon_names),
     target_names{it} = [W.gene_names{it} ' (' W.operon_names{it} ')'];
    end
end

blocks = zeros(size(W.data));

if options.blocks,
  dum = W.data;
  dum(:,find(sum(dum==0)==0)) = 0;
  dist = dum * dum' ;
  g.matrix = double(dist>0);
  [subgraphs,node_order] = find_separate_subgraphs(g);
  for it =1:length(subgraphs)
    blocks(subgraphs{it},find(sum(W.data(subgraphs{it},:),1)))=1;
  end
end
  
switch method,
  case 'M',  im(W.data,[],target_names,W.TF_names); 
  case 'R',  imagesc(W.data,[-1,1]); set(gca,'XTick',[],'YTick',[]); colormap(my_colors);
  case 'N',  
    set(gca,'FontSize',options.textsize);
    W.signs = W.signs +  0.7 * ( W.data .* ( W.signs==0) );
    if options.blocks, W.signs = W.signs + 0.3 * (blocks.*(W.data==0)); end
    imagesc(W.signs,[-1,1]); 
    set(gca,'FontSize',options.textsize,'YTick',1:size(W.data,1), 'YTicklabel',target_names); 
    if length(options.position),     set(gca,'Position',options.position); end
    my_xticklabel(1:length(W.TF_names),0.4,W.TF_names,2*options.textsize)
    colormap([red; red; red; 1 1 1; 1 1 1; 0.9 0.9 0.9; grey; blue]);
    
  otherwise, %% network representation

      if isfield(options,'W_support'),
        M_support = nca_network_support(W,options.W_support);
      else,
        M_support = zeros(size(W.data));
      end
    
      pos_operon  = length(target_names):-1:1;
      pos_TF_mean = sum(diag(length(target_names)+1:-1:2)*W.data,1) ./ (sum(W.data,1)+2);

      if length(W.TF_names)==1,
	pos_TF_distr = mean(pos_operon); 
      else,
	pos_TF_distr = 1 + [ 1-0.75/(length(W.TF_names)) : -1/(length(W.TF_names)) : 0.25/(length(W.TF_names))] * length(target_names);
%        if length(W.TF_names) < length(W.gene_names), pos_TF_distr = pos_TF_distr + 0.5 * pos_TF_distr(end-1); end
      end
      
      pos_TF = 0 .* pos_TF_mean + 1 * pos_TF_distr;
      
      [I0,J0]= ind2sub(size(W.data),find(W.data));      
      [nlinks,ord] = sort(-sum(W.data,1));     
      % arrange order such that TF with few connections are plotted last
      I = []; J = [];
      for it = 1:length(ord)
        ind = find(J0==ord(it));
        I = [I; I0(ind)];
        J = [J; J0(ind)];
      end

      if length(options.A_values),
        show_A = 1;
        cmap = rb_colors(250);
        if options.normalise_A,
          %% normalisation for entire A matrix; sqrt of max from both
          %% sides (same procedure as in W_display)          
          Amaxcol = nanmax(abs(options.A_values));
          Amaxcol(isnan(Amaxcol)) = 1;
          Amaxcol(Amaxcol==0) = 10^-10;
          Amaxrow = nanmax(abs(options.A_values'));
          Amaxrow(isnan(Amaxrow)) = 1;
          Amaxrow(Amaxrow==0) = 10^-10;
          options.A_values   = diag(1./sqrt(Amaxrow)) * options.A_values * diag(1./sqrt(Amaxcol));
          %% normalise A values to maximum absolute value per TF
          % Amax = nanmax(abs(options.A_values));
          % Amax(isnan(Amax)) = 1;
          % Amax(Amax==0) = 10^-10;
          % options.A_values = options.A_values * diag(1./Amax);
        elseif max(abs(options.A_values(:)))>1,
          error('A values out of range'); 
        end
        color_index    = ceil(1 + 248 * 0.5*[options.A_values+1]);
      end
      for it = 1:length(I),
        my_support = M_support(I(it),J(it));
	h =  line([options.xTF,options.xOP],[pos_TF(J(it)),pos_operon(I(it))],'LineWidth',options.lw + my_support);
        if show_A,
%          options.A_values(I(it),J(it))
          color = cmap(color_index(I(it),J(it)),:);
        else,
          switch W.signs(I(it),J(it)),
            case 1, color  = blue;
            case 0, color  = grey;
            case -1, color = red;
          end
        end
	set(h,'Color',color);

      %   if isfield(options,'A_quality'),
      %     if options.A_quality(I(it),J(it)),
      %       set(h,'LineStyle','--');
      %     end
      %   end
      
      end
      
      if ~options.no_gene_labels, 
      text(options.xOP*ones(size(pos_operon))+options.textOPxoffset,...
           pos_operon+options.textyoffset, target_names,'FontSize',options.textsize);
      end
      text(options.xTF*ones(size(pos_TF))+options.textTFxoffset,...
           pos_TF+options.textyoffset, W.TF_names,'FontSize',options.textsize);
      set(gca, 'XTick',[]); set(gca, 'YTick',[]);
      
      if ~isempty(options.arrow_values),
        for it1=1:size(W.data,1),
          for it2=1:size(W.data,2),
            if W.data(it1,it2)~=0,
              b= sum(W.data(it1,:));
              a= sum(W.data(:,it2));
              rndn = 0.4 + 0.3 * (a-b)/(a+b) + 0.1*randn;
              text( (options.xTF+0.1) + (rndn)*(1-options.xTF-0.1), ...
                    rndn*pos_operon(it1)+(1-rndn)*pos_TF(it2), num2str(options.arrow_values(it1,it2),2));
            end
          end
        end
      end
      
end

if show_A,
  set(gcf,'Color',[0.9 0.9 0.9],'InvertHardCopy', 'off');
end

%axis equal
axis tight
