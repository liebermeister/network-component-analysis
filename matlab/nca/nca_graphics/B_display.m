function B_display(B,TFlist,D,W,ni,color,t_values,user_pars,xlabellist)

eval(default('color','[0 0 0]','user_pars','struct','xlabellist','TFlist'));

pars.fill    = 1;
pars.textflag   = 0;
pars.print_flag = 0;
pars.fontsize   = [];
pars.linewidth  = 1;
pars.ytick = [];
pars.yline = [];
pars.print_experiments = 0;
pars.y_scale_only_once = 0;

pars = join_struct(pars,user_pars);
if isempty(pars.fontsize),
  switch pars.print_flag,
    case 0, pars.fontsize = 10;
    case 1, pars.fontsize = 4;
  end
end

if iscell(B),  
  B_std = B{2};  B = B{1}; 
else, 
  B_std = zeros(size(B));
end

if ~exist('D','var'), D = []; end
if isempty(D),
  D.experiments = {1:size(B,2)}; D.experiments = {''};
end

if ~exist('ni','var'),    ni = length(TFlist); end
if ~exist('color','var'), color = [ 0 0 1]; end
if ~exist('t_values','var'), t_values = 1:size(B,2); end

if length(ni)==1, ni = [ni,1]; end 

for it = 1:length(D.experiments),
  starting_points(it)= D.experiments{it}(1);
end

for it = 1:min(length(TFlist),prod(ni)),

  subplot(ni(1),ni(2),it);
  TFind = find_TF(TFlist(it),W); 
  hold on;

  if length(pars.yline),
    for it2 = 1:length(pars.yline),
      line( [0 t_values(end)], pars.yline(it2) * [1 1],'color',[0.5 0.5 0.5]);
    end
  end
  
  for it2 = 1:length(D.experiments),
    if length(D.experiments{it2}),
      if exist('B_std','var'),
        if pars.fill,
          plot_range(t_values(D.experiments{it2}),B(TFind,D.experiments{it2})',B_std(TFind,D.experiments{it2})',[],color);
        else,
          plot_range(t_values(D.experiments{it2}),B(TFind,D.experiments{it2})',B_std(TFind,D.experiments{it2})',[],color,color,color,0);
        end
      else,      
        plot(t_values(D.experiments{it2}),B(TFind,D.experiments{it2})','color',color,'LineWidth',pars.linewidth);
      end
    end
  end

  if length(pars.ytick),
    if [pars.y_scale_only_once==0] + [mod(it-1,ni(2))==0],
      set(gca,'YTick',pars.ytick);
    else,
      set(gca,'YTick',[]);
    end
    if length(pars.yline),
    mmin = min(min(pars.ytick), min(pars.yline));
    mmax = max(max(pars.ytick), max(pars.yline));
    else
    mmin = min(min(pars.ytick));
    mmax = max(max(pars.ytick));
    end
    axis([0 t_values(end) mmin mmax]);
    set(gca,'XTick',[]);
  else
    set(gca,'XTick', [],'YTick',[]);
    axis tight
    ax = axis; mmin = ax(3); mmax = ax(4);
  end

  for it2 = 1:length(D.experiments),
    line( t_values(starting_points(it2)) * [1 1], [mmin,mmax],'color',[0 0 0]);
  end

  hold off;
  if length(xlabellist),  title(xlabellist{it},'Fontsize',pars.fontsize); end

  if isfield(D,'ylabel'), if mod(it-1,ni(2))==0, ylabel(D.ylabel); end; end
  if isfield(D,'xlabel'), if it>[ni(1)-1]*ni(2), 
      xlabel(D.xlabel);       
      if pars.print_experiments, set(gca,'XTick',t_values(starting_points),'XTickLabel',D.experiment_names); end
    end; end
  
set(gca,'FontSize',pars.fontsize) 

end
