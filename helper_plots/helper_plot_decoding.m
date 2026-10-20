function helper_plot_decoding(cfg)
% this funciton plots decoding curves
% cfg can contain following fields:
%       * yavg: plain lines, dim1 = nsmp from epoch ,dim2 = nvar to plot
%       * yerr: shaded areas sem, dim1 = nsmp from epoch, dim2 = nvar to plot
%       * xvec: time samples
%       * rgb:  rgb color triplets for each variable to plot
%       * win:  gray-shaded area around peak
%       * peak: add peak vertical line

if isfield(cfg,'pbar'),pbar = cfg.pbar;
else
    pbar = 1.1;
end
if isfield(cfg,'HGT'),HGT = cfg.HGT;
else
    HGT = 4;
end

% decoded variables to plot in color, plain lines
if isfield(cfg,'yavg')
    yavg = cfg.yavg;
    nvar = size(yavg,2);
    nsmp = size(yavg,1);
else
    warning('nothing to plot!')
    return
end
if isfield(cfg,'rgb'),rgb = cfg.rgb;
else, rgb = lines; rgb = rgb(1:nvar,:);
end
if isfield(cfg,'yavg')&&(nvar>size(rgb,1))
    rgb = cat(1,rgb,lines);
end
if isfield(cfg,'yavg_dashed')
    if isfield(cfg,'rgb_dashed'), rgb_dashed = cfg.rgb_dashed;
    else, rgb_dashed = rgb;
    end
end

if isfield(cfg,'issubplot')
    if ~cfg.issubplot
        if isfield(cfg,'figname')
            figure('Color','white','Name',cfg.figname);
        else
            figure('Color','white');
        end
    end
else
    if isfield(cfg,'figname')
        figure('Color','white','Name',cfg.figname);
    else
        figure('Color','white');
    end
end

hold on

if isfield(cfg,'xvec'),xvec = cfg.xvec;
    nsmp = numel(xvec);
else
    xvec = 1:nsmp;
end

if isfield(cfg,'xlim'),xlim(cfg.xlim);
else
    xlim([-0.2 0.8]);
end

if isfield(cfg,'ylim'),yl = cfg.ylim;
else
    yl = [-0.025 0.4];
end
ylim(yl);

%% gray time window (averaging window)
if isfield(cfg,'win'), win = cfg.win;
    patch([min(win) max(win) max(win) min(win)],[yl(1) yl(1), yl(2) yl(2)],[.8 .8 .8],'EdgeColor','none','FaceAlpha',.5);
end

%% SHADED AREAS around lines
if isfield(cfg,'yerr')
    yerr = cfg.yerr;
    nvar_err = size(yerr,2);
    nvar = max(nvar,nvar_err);
    yerr_vec = nan(2,nsmp,nvar);
    for ivar = 1:nvar
        yerr_vec(1,:,ivar) = reshape(yavg(:,ivar),1,nsmp) + reshape(yerr(:,ivar),1,nsmp);
        yerr_vec(2,:,ivar) = reshape(yavg(:,ivar),1,nsmp) - reshape(yerr(:,ivar),1,nsmp);
        patch([xvec,fliplr(xvec)],[yerr_vec(1,:,ivar),fliplr(yerr_vec(2,:,ivar))], ...
            0.5*(rgb(ivar,:)+1),'EdgeColor','none','FaceAlpha',.5);
    end
end
if isfield(cfg,'yerr_dashed') && isfield(cfg,'yavg_dashed')
    yerr_dashed = cfg.yerr_dashed;
    yavg_dashed = cfg.yavg_dashed;
    nvar_err = size(yerr_dashed,2);
    yerr_vec = nan(2,nsmp,nvar_err);
    for ivar = 1:nvar_err
        yerr_vec(1,:,ivar) = reshape(yavg_dashed(:,ivar),1,nsmp) + reshape(yerr_dashed(:,ivar),1,nsmp);
        yerr_vec(2,:,ivar) = reshape(yavg_dashed(:,ivar),1,nsmp) - reshape(yerr_dashed(:,ivar),1,nsmp);
        patch([xvec,fliplr(xvec)],[yerr_vec(1,:,ivar),fliplr(yerr_vec(2,:,ivar))], ...
            0.5*(rgb_dashed(ivar,:)+1),'EdgeColor','none','FaceAlpha',.5);
    end
end
%% black lines
plot([0,0],ylim,'k:','LineWidth',0.75);
plot(xlim,[0,0],'Color',[.8 .8 .8],'LineWidth',0.75);

%% peak value
if isfield(cfg,'peak')
    xpeak = cfg.peak;
    plot(xpeak*[1 1],ylim,'k','LineWidth',0.5);
end
    
%% decoded variables
for ivar = 1:nvar
    plot(xvec,yavg(:,ivar),'Color',rgb(ivar,:),'LineWidth',1.5);
    if isfield(cfg,'peak')
        xpeak = cfg.peak;
        % add on time course
        ypeak = interp1(xvec,yavg(:,ivar),xpeak);
        plot(xpeak,ypeak,'o','MarkerFaceColor',rgb(ivar,:),'MarkerSize',3,'MarkerEdgeColor','w');
    end
end
if isfield(cfg,'yavg_dashed')
    yavg_dashed = cfg.yavg_dashed;
    nvar_dashed = size(yavg_dashed,2);
    for ivar = 1:nvar_dashed
        plot(xvec,yavg_dashed(:,ivar),'Color',rgb_dashed(ivar,:),'LineWidth',1.5,'LineStyle',':');
    end
end

%% stats
if isfield(cfg,'stats')
    if isfield(cfg.stats,'h')
        h = cfg.stats.h;
        nstats = size(h,1);
        if isfield(cfg.stats,'rgb')
            rgb_stats = cfg.stats.rgb;
        else
            rgb_stats = rgb;
        end
        if size(rgb_stats,1)<nstats
            rgb_stats = cat(1,rgb_stats,lines);
        end
        for istats = 1:nstats
            c = bwconncomp(h(istats,:)==1);
            for ic = 1:c.NumObjects
                ifirst = c.PixelIdxList{ic}(1);
                ilast  = c.PixelIdxList{ic}(end);
                plot([xvec(ifirst) xvec(ilast)],[max(ylim) max(ylim)],'Color',rgb_stats(istats,:),'Linewidth',2);
            end
        end
    end
end

%%

if isfield(cfg,'xtick'), xt = cfg.xtick;
else
    xt = -.2:.2:.8;
end
if isfield(cfg,'xticklabel'), xtl = cfg.xticklabel;
else
    xtl = {'-0.2' '0' '0.2' '0.4' '0.6' '0.8'};
end
if isfield(cfg,'ytick')
    set(gca,'YTick',cfg.ytick);
end

if isfield(cfg,'xlabel'), xlab = cfg.xlabel;
else
    xlab = 'time from stimulus onset (s)';
end
if isfield(cfg,'ylabel'), ylab = cfg.ylabel;
else
    ylab = 'coding precision (z)';
end

set(gca,'XTick',xt,'XTickLabel',xtl);
xlabel(xlab)
ylabel(ylab)
helper_ILL(gcf,pbar,HGT);

if isfield(cfg,'saveornot') && isfield(cfg,'figname') && isfield(cfg,'figpath')
    if cfg.saveornot
        print(fullfile(cfg.figpath,cfg.figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',cfg.figname,cfg.figpath)
    end
end
        