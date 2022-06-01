function hfig = plot_surface_movie(p,sourcedata,parameter,surface_inflation,single_plot,interptype)
    % Render volume data onto a surface
    %
    % *REQUIRES WORKBENCH*
    %
    % INPUTS
    % - p - parcellation object
    % - data - data that can be saved to a nii file via p.savenii(data)
    % - surface_inflation - integer level of inflation for display surface
    %   (default=0, no inflation)
    % - single_plot - put the two hemispheres together (default=false, like
    %   Workbench)
    % - interptype - passed to workbench (default='trilinear')
    %
    % If the data has more than one volume, the output plot will also have
    % more than one volume. The figure has a property 'current_vol' that lets
    % you set the which one is displayed e.g.
    %
    % h = p.plot_surface(data) h.current_vol = 2; % Display second volume
    %
    % The buttons on the top of the figure let you manually cycle through
    % volumes, or display them in sequence
    sz = size(sourcedata.(parameter));
    data=reshape(sourcedata.(parameter), sz(1), []);

    if nargin < 6 || isempty(interptype) 
        interptype = 'trilinear';
    end
    
    if nargin < 5 || isempty(single_plot) 
        single_plot = true;
    end

    if nargin < 4 || isempty(surface_inflation) 
        surface_inflation = 0;
    end

    if single_plot == true && surface_inflation ~= 0
        fprintf(2,'Warning, single plot with inflated surface does not currently plot correctly');
    end
    
    niifile = p.savenii(data);
    output_right    = [niifile '_right.func.gii'];
    output_left     = [niifile '_left.func.gii'];
    cl = onCleanup(@()  cellfun(@delete,{niifile,output_left,output_right})); % Enable deleting temp files even if debugging

    surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
    surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');

    switch surface_inflation
        case 0
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');
        case 1
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.inflated.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.inflated.32k_fs_LR.surf.gii');
        case 2
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.very_inflated.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.very_inflated.32k_fs_LR.surf.gii');
    end


    % Map volume to surface
    runcmd('wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_right,output_right,interptype)
    runcmd('wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_left,output_left,interptype)

    sl = gifti(display_surf_left);
    vl = gifti(output_left);
    sr = gifti(display_surf_right);
    vr = gifti(output_right);


brainordinate=[];
if single_plot
  brainordinate.pos=[sl.vertices; sr.vertices];
  
else
  brainordinate.pos=[sl.vertices; [1 -1 1].*(sl.vertices + 1.3*[0,range(sl.vertices(:,2)),0])];
end
brainordinate.tri=[sl.faces; max(sl.faces)+sr.faces];
for k=1:38
  brainordinate.brodmannlabel{k+1} = sprintf('ROI %d', k);
end
brainordinate.brodmannlabel{1} = 'None';
sourcedata.label=brainordinate.brodmannlabel;
load('Giles39_atlas.mat')
brainordinate.brodmann = atlas.brodmann;
sourcedata.pos=brainordinate.pos;
sourcedata.tri=brainordinate.tri;
sourcedata.(parameter) = [vl.cdata; vr.cdata];
cfgp.funparameter=parameter;
if numel(size(sourcedata.(parameter)))<=2
  cfgp.xparam='freq';
  cfgp.yparam=[];
else
  cfgp.xparam='time';
  cfgp.yparam='freq';
end

cfgp.funcolormap = flipud(brewermap(16, 'RdBu'));
% cfgp.funcolorlim = [-2.5 2.5];
% cfgp.colormap=cfgp.funcolormap;
ft_sourcemovie(cfgp, sourcedata)