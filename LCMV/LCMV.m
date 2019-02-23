%LCMV method from Leadfield
%http://www.fieldtriptoolbox.org/tutorial/aarhus/beamformingerf?s[]=lcmv#compute_leadfield_at_location_m1
%Qinyuan Wei 2017.8
%% load models and templetes
load data_eeg_reref_ica.mat
load headmodel_eeg.mat
load elec.mat
load mri_segmented.mat
%% sort into left and right hand response
cfg.trials       = find(data_eeg_reref_ica.trialinfo(:,1) == 4096);
data_eeg_right    = ft_redefinetrial(cfg, data_eeg_reref_ica);
%%
cfg = [];
cfg.toilim = [-.15 -.05];
datapre = ft_redefinetrial(cfg, data_eeg_right);
cfg.toilim = [.05 .15];
datapost = ft_redefinetrial(cfg, data_eeg_right);
%% output cov matrix of the entire interval
cfg = [];
cfg.covariance='yes';
cfg.covariancewindow = [-.15 .15];
cfg.vartrllength = 2;
avg = ft_timelockanalysis(cfg,data_eeg_right);
cfg = [];
cfg.covariance='yes';
avgpre = ft_timelockanalysis(cfg,datapre);
avgpst = ft_timelockanalysis(cfg,datapost);

cfg                 = [];
cfg.elec         = elec;
cfg.channel         = data_eeg_right.label;
cfg.vol             = headmodel_eeg;
cfg.dics.reducerank = 3; % default for MEG is 2, for EEG is 3
cfg.grid.resolution = 8;   % use a 3-D grid with a 1 cm resolution, the resolution can be changed
cfg.grid.unit       = 'mm';%unit can be changed
cfg.grid.tight      = 'yes';
cfg.normalize = 'yes';
grid = ft_prepare_leadfield(cfg);
save grid_eeg grid

cfg=[];
cfg.method='lcmv';
cfg.grid=grid;%leadfield
cfg.elec = elec;
cfg.vol=headmodel_eeg;
cfg.lcmv.keepfilter='yes';
cfg.lcmv.lambda = '5%';
cfg.channel = data_eeg_right.label;
cfg.senstype = 'EEG';
sourceavg=ft_sourceanalysis(cfg, avg);

%%
cfg=[];
cfg.method='lcmv';
cfg.elec = elec;
cfg.grid=grid;
cfg.grid.filter=sourceavg.avg.filter;
cfg.vol=headmodel_eeg;
cfg.lcmv.lambda = '5%';
cfg.channel           = data_eeg_right.label;
cfg.senstype = 'EEG';
sourcepreM1=ft_sourceanalysis(cfg, avgpre);
sourcepstM1=ft_sourceanalysis(cfg, avgpst);

M1eeg=sourcepstM1;
M1eeg.avg.pow=(sourcepstM1.avg.pow-sourcepreM1.avg.pow)./sourcepreM1.avg.pow;

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, M1eeg, mri_segmented);
%%
source_int.mask = source_int.pow > max(source_int.pow(:))*.3;% a threshod which decides the present of sources
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.location = [-28 -17 67]; 
cfg.funcolorlim = [-.2 .2];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg,source_int);

figure(2);
title('Result');
ft_plot_vol(headmodel_eeg, 'edgecolor', 'none', 'facealpha', 0.2);
hold on;
gridin=M1eeg.pos(M1eeg.inside,:);
scatter3(gridin(:,1),gridin(:,2),gridin(:,3),'filled','cdata',abs(M1eeg.avg.pow(M1eeg.inside)));
colorbar;
colormap('jet')

figure(3);
title('Sourceavg');
ft_plot_vol(headmodel_eeg, 'edgecolor', 'none', 'facealpha', 0.2);
hold on;
gridin=sourceavg.pos(sourceavg.inside,:);
scatter3(gridin(:,1),gridin(:,2),gridin(:,3),'filled','cdata',abs(sourceavg.avg.pow(sourceavg.inside)));
colorbar;
colormap('jet')
