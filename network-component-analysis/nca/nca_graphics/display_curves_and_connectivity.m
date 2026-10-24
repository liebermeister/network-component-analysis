function options = display_curves_and_connectivity(W,X,TF_act,X_additional,options,arrow_values);

% display_curves_and_connectivity(W,X,TF_act,X_additional,options,arrow_values);
%options_default = struct('arrow_numbers',1,'implot',0,'lw',1,'net_ratio',0.3,'box_ratio', 0.25,'textyoffset',-0.25);

options_default = struct('arrow_numbers',1,'implot',0,'lw',1,'net_ratio',5, 'box_ratio', 0.25,'textyoffset',0,'linepositions',[],'position',[],'net_width',[],'individual_scaling',1,'no_gene_labels',0);

colors = BEcolor;

n_genes = length(W.gene_names);
n_TF    = length(W.TF_names);

if ~exist('options','var'), options = []; end
if isempty(options), options = struct; end

ff = fields(options_default);
for it = 1:length(ff),
  if ~isfield(options,ff{it}),
    options=setfield(options,ff{it},getfield(options_default,ff{it}));
  end
end

if ~exist('TF_act','var'),             TF_act = []; end
if ~exist('X_additional','var'), X_additional = []; end
if ~exist('arrow_values','var'), arrow_values = []; options.arrow_numbers =0; end

if isempty(options.net_width),  options.net_width = sqrt(1+length(W.TF_names)); end

if ~isfield(options,'textsize'),
options_t         = ceil(18/[sqrt(length(W.gene_names))]);
options.textsize  = options_t;
end

%if isempty(options.net_width), options.net_width = size(X,1) / options.net_ratio; end 
options.boxwidth = 1 / options.box_ratio;

if ~isfield(options,'xTF'), options.xTF = 0;                   end
if length(options.net_width),
  options.xOP = options.xTF + options.net_width;
end

if ~isfield(options,'label_width'),
  options.label_width = max(0.3 * options.net_width, 0.2 * options.textsize);
end

if ~isfield(options,'textTFxoffset'), options.textTFxoffset = -1.5  * options.label_width; end 
if ~isfield(options,'textOPxoffset'), options.textOPxoffset = 0.5   * options.label_width; end 
if ~isfield(options,'boxTFxoffset'),  options.boxTFxoffset  = options.xTF - 2 * options.label_width - options.boxwidth; end 
if ~isfield(options,'boxOPxoffset'),  options.boxOPxoffset  = options.xOP + 2 * options.label_width; end 
if options.no_gene_labels, options.boxOPxoffset  = options.xOP + 0.2 * options.label_width; end
  
if ~isfield(options,'time_values'),   options.time_values   = 1:size(TF_act,2); end

% normalise options.time_values:

options.time_values = options.time_values - options.time_values(1);
options.time_values = options.time_values / options.time_values(end);

% --------------------------------------------------------------------
% draw network

clf;
[pos_TF,pos_operon] = W_display(W,'network',options); 
hold on;

if ~isempty(TF_act),
  if iscell(TF_act), 
    TF_act_rest = TF_act(2:end);
    TF_act      = TF_act{1};
  else, 
    TF_act_rest = [];
  end
end

% --------------------------------------------------------------------
% TF curves

xvalues   = options.boxTFxoffset + options.boxwidth * options.time_values;
xline     = options.boxTFxoffset + options.boxwidth * options.linepositions/(size(TF_act,2)-1);
fullysize = 0.8*(2+size(W.data,1))/(2+size(W.data,2));

if ~isempty(TF_act),

  hold on;
  
  if options.implot,
    h = imagesc(TF_act,max(max(abs(TF_act)))*[-1 1]); 
    set(h,'XData',[xvalues(1), xvalues(end)],'YData',[pos_TF(1),pos_TF(end)]); colormap(rb_colors);
    hold off;
  else
    
    for it = 1:length(options.linepositions), 
      line([xline(it),xline(it)],[0.5,n_genes+.5],'color',[.7 .7 .7],'LineWidth',options.lw); hold on;
    end
    
    if options.individual_scaling,
      %% normalisation per TF
      %%      TF_act = TF_act - repmat(mean(TF_act,2),1,size(TF_act,2));
      %%      if length(TF_act_rest), % range
      %%        TF_act_rest{1} = TF_act_rest{1}-repmat(mean(TF_act_rest{1},2),1,size(TF_act_rest{1},2));
      %%        TF_act_rest{2} = TF_act_rest{2}-repmat(mean(TF_act_rest{2},2),1,size(TF_act_rest{2},2));      
      %%  normalisation for all data
      maxv = max(max(abs(TF_act')),0.05)';
      %%      maxv = max(abs(TF_act'))';
     if length(pos_TF)>1,
      delta_posTF = abs(pos_TF(2) - pos_TF(1));
     else,
      delta_posTF = 1;
     end
     ysize = [0.4 ./ maxv] * delta_posTF * ones(1,size(TF_act,1));
    else 
     ysize = ones(1,size(TF_act,1));
    end

    if length(TF_act_rest),
    for it = 1:n_TF,
      vv = [ pos_TF(it) + ysize(it)*( TF_act_rest{1}(it,:)); ...
             pos_TF(it) + ysize(it)*( TF_act_rest{2}(it,:))];
      simple_range_plot( xvalues', vv, colors.TFlight);
    end
    end
    
    for it = 1:n_TF,
      plot( xvalues, pos_TF(it) + ysize(it)*(TF_act(it,:)),'Color',colors.TF,'LineWidth',options.lw);
    end
    if length(TF_act_rest)>2,
      for it2 = 3:length(TF_act_rest),
        for it = 1:n_TF,
          plot( xvalues, pos_TF(it) + ysize(it)*(TF_act_rest{it2}(it,:)),'LineWidth',options.lw,'color',[0 0 0]);
        end
      end
    end
    
  end

end

% --------------------------------------------------------------------
% gene curves

xvalues = options.boxOPxoffset +  options.boxwidth * options.time_values;
xline   = options.boxOPxoffset +  options.boxwidth * ( options.linepositions/(size(TF_act,2)-1));

hold on;  

if options.implot,
  h = imagesc(flipud(X),max(max(abs(X)))*[-1 1]); set(h,'XData',[xvalues(1) xvalues(end)]); colormap(rb_colors);
  hold off;
else,
  
  for it = 1:length(options.linepositions),
    line([xline(it),xline(it)],[0.5,n_genes+.5],'color',[.7 .7 .7],'LineWidth',options.lw); hold on;
  end
  
  if options.individual_scaling,
    %% normalisation per gene
    if size(X_additional),
      if iscell(X_additional),
        %% assume that the first entry is a mean and the others are
        %lower and upper bounds
        for it = 1:length(X_additional),
          X_additional{it} = X_additional{it}-repmat(mean(X,2),1,size(X,2));
        end
      else,
        X_additional = X_additional - repmat(mean(X,2),1,size(X,2));
      end
    end
    X = X-repmat(mean(X,2),1,size(X,2));
    maxv  = max(max(X') - min(X'),1);
    ysize = 0.6./maxv;
  else,
    %%  normalisation for all data
    maxv = max(max(abs(X)));
    ysize = 0.7/maxv*ones(1,size(X,1));
  end
  
  for it = 1:n_genes,
   if ~isempty(X_additional), 
     if ~iscell(X_additional),
       X_additional = {X_additional};
     end
   end 

   if ~isempty(X_additional),
     if length(X_additional)>1,
       vv = [ pos_operon(it) + ysize(it)*( X_additional{2}(it,:)); ...
              pos_operon(it) + ysize(it)*( X_additional{3}(it,:))];
       simple_range_plot( xvalues', vv, colors.transcript_light);
     end
     plot(xvalues, pos_operon(it) + ysize(it)*(X_additional{1}(it,:)),'LineWidth',options.lw,'color',colors.transcript_light); hold on;
   end
  
   %% actual values
   plot(xvalues, pos_operon(it) + ysize(it)*(X(it,:)),'Color',colors.transcript,'LineWidth',options.lw); hold on;

  end
end


hold off;
axis equal
axis tight
axis off

function simple_range_plot(t,x,col)

if size(t,2) ~= 1, t = t'; end
if size(x,2) ~= 2, x = x'; end

h = fill([t;flipud(t)],[x(:,1); flipud(x(:,2))],col,'EdgeColor',col);
