function helper_plot_source(cfg)
% map on the cortex of the default anatomy, one figure per view: replaces the display from the
% Brainstorm database (script_view_sources) in the figure scripts, with the display rules of Brainstorm
% (figure_3d.m: ThresholdSurfaceData, BlendAnatomyData, lights and surface properties; default display
% options: amplitude threshold 50%, sulci shown; smoothing 40% as in the figures of the paper)
%   cfg.cortex     cortex of the default anatomy (default_cortex.mat: Vertices, Faces, SulciMap)
%   cfg.map        template vertex x 1, values to display
%   cfg.clim       color limits [min max] (SetMaxCustom in Brainstorm)
%   cfg.cmap       colormap (n x 3)    [default: Brainstorm colormap for the sources, cmap_royal_gramma]
%   cfg.absolute   absolute values (SetColormapAbsolute)  [default: false]
%   cfg.threshold  amplitude threshold (fraction of the color range, Brainstorm "Amplitude" slider): vertices
%                  below it are not colored  [default: 0.5, Brainstorm default]
%   cfg.views      cell array of 'top', 'bottom', 'left', 'right', 'front', 'back'  [default: {'top','left'}]
%   cfg.smooth     smoothing of the cortex for the display (Brainstorm "Smooth" slider, 0 to 1, e.g. 0.3 for 30 percent,
%                  same computation as figure_3d.m: tess_smooth with ceil(300*smooth*nvertices/1e5) iterations)
%                  [default: 0.4, as in the figures of the paper made with Brainstorm]
%   cfg.figname    figure name (the view is added: figname_top.pdf)
%   cfg.figpath, cfg.saveornot

if ~isfield(cfg,'cmap')
    cfg_path = config_paths;
    addpath(fullfile(cfg_path.brainstorm,'toolbox','misc')) % Brainstorm colormaps
    cfg.cmap = cmap_royal_gramma(256);
end
if ~isfield(cfg,'absolute'), cfg.absolute = false; end
if ~isfield(cfg,'threshold'), cfg.threshold = 0.5; end
if ~isfield(cfg,'views'), cfg.views = {'top','left'}; end
if ~isfield(cfg,'smooth'), cfg.smooth = 0.4; end
if ~isfield(cfg,'saveornot'), cfg.saveornot = false; end

V = cfg.cortex.Vertices;
F = cfg.cortex.Faces;
% smoothing of the cortex for the display (SmoothSurface in figure_3d.m of Brainstorm)
if cfg.smooth > 0
    cfg_path = config_paths;
    addpath(genpath(fullfile(cfg_path.brainstorm,'toolbox'))) % tess_vertconn, tess_smooth
    niter = ceil(300*cfg.smooth*size(V,1)/100000);
    V = tess_smooth(V,cfg.smooth,niter,tess_vertconn(V,F),1);
end
data = double(cfg.map(:));
data(isnan(data)) = 0;
if cfg.absolute
    data = abs(data);
end
lim = cfg.clim;

% amplitude threshold (ThresholdSurfaceData): thresholded vertices set to 0, i.e. not colored
if ~cfg.absolute && lim(1) == -lim(2)
    data(abs(data) < cfg.threshold*max(abs(lim))) = 0;
elseif lim(2) <= 0
    data(data < lim(1)+(lim(2)-lim(1))*cfg.threshold) = lim(1);
    data(data > lim(2)) = 0;
else
    data(data < lim(1)+(lim(2)-lim(1))*cfg.threshold) = 0;
    data(data > lim(2)) = lim(2);
end

% colors (BlendAnatomyData): cortex in grey (gyri 0.6, sulci 0.45), colormap on the non-zero vertices
anatcol = [.45 .45 .45; .6 .6 .6];
if isfield(cfg.cortex,'SulciMap')
    C = anatcol(2-double(cfg.cortex.SulciMap(:)),:);
else
    C = repmat(anatcol(2,:),size(V,1),1);
end
ncol = size(cfg.cmap,1);
icol = round((ncol-1)/(lim(2)-lim(1))*(data-lim(1)))+1;
icol = min(max(icol,1),ncol);
iblend = find(data ~= 0);
C(iblend,:) = cfg.cmap(icol(iblend),:);

% views (SetStandardView): camera direction and up vector in the coordinates of the anatomy
camdir = struct('top',[0 0 1],'bottom',[0 0 -1],'left',[0 1 0],'right',[0 -1 0],'front',[1 0 0],'back',[-1 0 0]);
camupv = struct('top',[1 0 0],'bottom',[1 0 0],'left',[0 0 1],'right',[0 0 1],'front',[0 0 1],'back',[0 0 1]);
% camera fixed for each view (independent of the size of the window on the screen): centre of the cortex,
% distance and zoom angle such that the whole cortex fits in the figure whatever the view
ctr = (min(V)+max(V))/2;
rad = max(sqrt(sum((V-ctr).^2,2))); % radius of the sphere containing the cortex
dst = 10*rad;
vang = 2*atand(1.02*rad/dst);
setview = @(vnam)set(gca,'CameraTarget',ctr,'CameraPosition',ctr+dst*camdir.(vnam), ...
    'CameraUpVector',camupv.(vnam),'CameraViewAngle',vang);
% one figure per map, the cortex drawn once, then one view after the other (SetStandardView in Brainstorm:
% camera and headlight moved); saved as in the previous figure scripts: figure of 180 x 180 pixels,
% print -vector -dpdf
figure('Color','white','Name',cfg.figname);
axes('Position',[0 0 1 1]); % cortex filling the figure, as in the 3D figures of Brainstorm
patch('Faces',F,'Vertices',V,'FaceVertexCData',C,'FaceColor','interp','EdgeColor','none', ...
    'BackFaceLighting','lit','AmbientStrength',0.5,'DiffuseStrength',0.5,'SpecularStrength',0.2, ...
    'SpecularExponent',1,'SpecularColorReflectance',0.5,'FaceLighting','gouraud');
axis equal off
colormap(gca,cfg.cmap);
set(gca,'CLim',lim); % for a colorbar
% lights of Brainstorm (figure_3d.m): 5 fixed lights placed from the top view, 1 headlight, intensity 0.4
setview('top');
hl = [camlight(0,40,'infinite'),camlight(180,40,'infinite'),camlight(0,-90,'infinite'), ...
    camlight(90,0,'infinite'),camlight(-90,0,'infinite'),camlight('headlight')];
set(hl,'Color',.4*[1 1 1]);
for iview = 1:numel(cfg.views)
    vnam = cfg.views{iview};
    figname = sprintf('%s_%s',cfg.figname,vnam);
    set(gcf,'Name',figname);
    setview(vnam);
    camlight(hl(6),'headlight'); % headlight follows the camera
    set(gcf,'Position',[0 0 180 180])
    drawnow
    if cfg.saveornot
        print(fullfile(cfg.figpath,figname),'-vector','-dpdf');
        fprintf('\nfigure %s saved in %s !\n',figname,cfg.figpath)
    end
end

end
