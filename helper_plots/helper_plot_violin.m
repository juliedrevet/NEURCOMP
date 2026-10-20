function helper_plot_violin(cfg)
if isfield(cfg,'HGT'),HGT = cfg.HGT;
else, HGT = 4 ;end
if isfield(cfg,'rgb'),rgb = cfg.rgb;
else, rgb = lines; end
if isfield(cfg,'rgb_dat'),rgb = cfg.rgb_dat;
else, rgb_dat = [0.5 0.5 0.5]; end
if isfield(cfg,'width'),wdth = cfg.width;
else, wdth = 0.3; end
if isfield(cfg,'pbar'),pbar = cfg.pbar;
else, pbar = 1 ;end
if isfield(cfg,'ylim'),yl = cfg.ylim;
else, yl = [0 1];end
if ~isfield(cfg,'ttest'),cfg.ttest = false;end
if ~isfield(cfg,'type'),cfg.type = 'vertical';end
if isfield(cfg,'xoff'),xoff = cfg.xoff; else, xoff = 0.6;end
if isfield(cfg,'ViolinEdgeColor'), ViolinEdgeColor = cfg.ViolinEdgeColor;
else ViolinEdgeColor = rgb;end
if isfield(cfg,'ViolinFaceColor'), ViolinFaceColor = cfg.ViolinFaceColor;
else ViolinFaceColor = rgb;end
if isfield(cfg,'plotmed'),plotmed = cfg.plotmed; else, plotmed = true;end
if isfield(cfg,'plotavg'),plotavg = cfg.plotavg; else, plotavg = false;end
if isfield(cfg,'plotsem'),plotsem = cfg.plotsem; else, plotsem = false;end
if isfield(cfg,'mrk_siz'),mrk_siz = cfg.mrk_siz; else, mrk_siz = 3;end
if isfield(cfg,'MarkerEdgeColor'), MarkerEdgeColor = cfg.MarkerEdgeColor;
else,MarkerEdgeColor = 'w';end


switch cfg.type
    case 'vertical'
        ori = 'v';
    case 'horizontal'
        ori = 'h';
end

if isfield(cfg,'x')&&(size(cfg.x,2)>size(rgb,1))
    rgb = cat(1,rgb,lines);
end

if isfield(cfg,'x')
    nmodel = size(cfg.x,2);
else
    return;
end

if isfield(cfg,'figname')
    figure('Color','white','Name',cfg.figname);
else    
    figure('Color','white');
end

if isfield(cfg,'xlabel'), xlb = cfg.xlabel;
else, xlb = 'Xlabel';end
if isfield(cfg,'ylabel'), ylb = cfg.ylabel;
else, ylb = 'Ylabel';end

if isfield(cfg,'labels'),labels = cfg.labels;
else, labels = cellfun(@(s)sprintf('Label%s',s),cellstr(num2str((1:nmodel)')),'UniformOutput',false);end

if isfield(cfg,'xticklabel'), xtlb = cfg.xticklabel;
else, xtlb = labels;end

if isfield(cfg,'xpos'), xpos = cfg.xpos;
else, xpos = 1:nmodel;end

hold on

if isfield(cfg,'x')
    x = cfg.x;
    if ~isfield(cfg,'ylim')
        if iscell(x)
            yl = [min(cellfun(@min,x)) max(cellfun(@max,x))];
        else
            yl = [min(x(:)) max(x(:))];
        end
    end
    if isfield(cfg,'pbar'),pbar = cfg.pbar;
    else, pbar = (nmodel+1.2)/HGT ;end
    if isfield(cfg,'xlim')
        xl = cfg.xlim;
    else
        xl = [1- xoff,nmodel+xoff];
    end
    switch ori
        case 'v'
            xlim(xl);
            ylim(yl);
        case 'h'
            ylim(xl);
            xlim(yl);
    end
    for imod = 1:nmodel
        if iscell(x)
            xdat = x{imod};
        else
            xdat = x(:,imod);
        end
        pos = xpos(imod);
        se = std(xdat)/sqrt(numel(xdat));
        xk = linspace(min(xdat),max(xdat),100);
        [pk,~,u] = ksdensity(xdat,xk);
        while true
            z = findpeaks(pk,xk);
            if numel(z) <= 1
                break
            else
                u = u*1.1;
                pk = ksdensity(xdat,xk,'BandWidth',u);
            end
        end
        pk = pk/max(pk);
        str = interp1(xk,pk,xdat);
        xmed = interp1(xk,pk,median(xdat));
        jit = linspace(-1,+1,numel(xdat));
        
        switch ori
            case 'v'
                if ~strcmp(ViolinFaceColor,'None')
                    patch([pos+pk*wdth,fliplr(pos-pk*wdth)],[xk,fliplr(xk)],0.5*(rgb(imod,:)+1),'EdgeColor','none','FaceAlpha',0.5);
                end
                if ~strcmp(ViolinEdgeColor,'None')
                    plot(pos+pk*wdth,xk,'-','Color',rgb(imod,:),'LineWidth',0.5);
                    plot(pos-pk*wdth,xk,'-','Color',rgb(imod,:),'LineWidth',0.5);
                end
                for i = 1:numel(xdat)
                    plot(xpos(imod)+jit(i)*str(i)*wdth,xdat(i),'o','MarkerFaceColor',.5*(rgb(imod,:)+1),'MarkerSize',mrk_siz,'MarkerEdgeColor',MarkerEdgeColor);
                end
                if plotsem
                    plot(pos*[1,1],mean(xdat)+se*[-1,+1],'k-','LineWidth',0.75);
                else
                    plot(pos*[1,1],mean(xdat)+std(xdat)*[-1,+1],'k-','LineWidth',0.75);
                end
                plot(pos,mean(xdat),'ko','MarkerFaceColor',rgb(imod,:),'MarkerSize',5,'LineWidth',0.75);
                if plotmed
                    plot(pos+xmed*wdth*[-1;+1],median(xdat).*[1 1],'Color',rgb(imod,:),'LineWidth',0.75);
                end
            case 'h'
                patch([xk,fliplr(xk)],[pos+pk*wdth,fliplr(pos-pk*wdth)],0.5*(rgb(imod,:)+1),'EdgeColor','none','FaceAlpha',0.5);
                plot(xk,pos+pk*wdth,'-','Color',rgb(imod,:),'LineWidth',0.5);
                plot(xk,pos-pk*wdth,'-','Color',rgb(imod,:),'LineWidth',0.5);
                for i = 1:numel(xdat)
                    plot(xdat(i),xpos(imod)+jit(i)*str(i)*wdth,'o','MarkerFaceColor',.5*(rgb(imod,:)+1),'MarkerSize',3,'MarkerEdgeColor','w');
                end
                plot(mean(xdat)+std(xdat)*[-1,+1],pos*[1,1],'k-');
                plot(mean(xdat),pos,'ko','MarkerFaceColor','k','MarkerSize',6,'LineWidth',0.5);
                plot(median(xdat).*[1 1],pos+xmed*wdth*[-1;+1],'Color',rgb(imod,:),'LineWidth',0.5);
        end

    end
    
    switch ori
        case 'v'
            set(gca,'XTick',xpos,'XTickLabel',xtlb(1:nmodel))
        case 'h'
            set(gca,'YTick',xpos,'YTickLabel',xtlb(1:nmodel))
    end
    
end


if isfield(cfg,'x') && cfg.ttest
    x = cfg.x;
    nmodel = size(x,2);
    if nmodel == 2
        if iscell(cfg.x)
            [H,P] = ttest2(x{1},x{2});
        else
            [H,P] = ttest(x(:,1),x(:,2));
        end
        switch ori
            case 'v'
                line([1 2],.9*max(ylim)*[1 1],'Color','k')
                if H == 0
                    text(1.5-.1,min(ylim)+.95*diff(ylim),'n.s.')
                else
                    if P<=.001
                        text(1.5-.1,min(ylim)+.95*diff(ylim),'***')
                    elseif P<=.01
                        text(1.5-.05,min(ylim)+.95*diff(ylim),'**')
                    elseif P<=.05
                        text(1.5-.01,min(ylim)+.95*diff(ylim),'*')
                    end
                    
                end
            case 'h'
                line(.9*max(xlim)*[1 1],[1 2],'Color','k')
                if H == 0
                    text(min(xlim)+.95*diff(xlim),1.5-.1,'n.s.')
                else
                    if P<=.001
                        text(min(xlim)+.95*diff(xlim),1.5-.1,'***')
                    elseif P<=.01
                        text(min(xlim)+.95*diff(xlim),1.5-.05,'**')
                    elseif P<=.05
                        text(min(xlim)+.95*diff(xlim),1.5-.01,'*')
                    end
                    
                end
        end
    end
end


if isfield(cfg,'ytick')
    switch ori
        case 'v'
            set(gca,'YTick',cfg.ytick);
    end
end
if isfield(cfg,'yticklabel')
    switch ori
        case 'v'
            set(gca,'YTickLabel',cfg.yticklabel);
        case 'h'
            set(gca,'XTickLabel',cfg.yticklabel);
    end
end
switch ori
    case 'v'
        xlim(xl);
        ylim(yl);
    case 'h'
        ylim(xl);
        xlim(yl);
end
switch ori
    case 'v'
        xlabel(xlb); ylabel(ylb);
    case 'h'
        ylabel(xlb); xlabel(ylb);
end

fig = helper_ILL(gcf,pbar,HGT,false);
set(fig,'PaperPositionMode','manual', ...
    'PaperPosition',[2.5,13,16,HGT],'PaperUnits','centimeters', ...
    'PaperType','A4','PaperOrientation','portrait');
if isfield(cfg,'saveornot') && isfield(cfg,'figname') && isfield(cfg,'figpath')
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end
    
end