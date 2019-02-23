%%
% MCE method in 3D grid model
% Qinyuan Wei 2017.8

%% load models and templetes
load('source_10240.mat');%ori
load ('standard_sourcemodel3d8mm.mat');%3D sourcemodel
load ('headmodel_template.mat');%vol
load ('elec_template.mat');%elec
%% compute leadfield
vol=ft_convert_units(vol, 'cm');%convert vol and sourcemodel into the same unit
sourcemodelwei = ft_convert_units(sourcemodel, 'cm');
cfg = [];
cfg.grid.pos = sourcemodelwei.pos;              % source points
cfg.grid.inside = 1:size(sourcemodelwei.pos,1); % all source points are inside of the brain
cfg.elec = elec;                      % sensor positions
cfg.vol = vol;                               % volume conduction model
cfg.normalize = 'yes';
cfg.normalizeparam  = 0.5;
K = ft_prepare_leadfield(cfg);
save K; 
%%
lf = K.leadfield;
num_channel = size(lf{1},1);
num_voxel = size(lf,2);
LFM=cell2mat(K.leadfield);
LFM_nrm = zeros(num_channel, num_voxel);
for i = 1:num_voxel
    LFM_nrm(:,i) = sqrt(LFM(:,(3*i-2)).^2+LFM(:,(3*i-1)).^2+LFM(:,(3*i)).^2);
end
%% generate ori
en = 2;         % electronic+electrode noise, 2~5uV
snr = 20;       % signal-to-noise-ratio

K = LFM_nrm; %K is the lead field matrix of the inverse problem

randn('seed',200); %the next randn would use the same series
ori=zeros(num_voxel,1);
ori(3283:3285,1)=0.5;
ori(3284,1)=1;
ori(2200,1)=1;
ori(2195,1)=0.5;
ori(2376,1)=0.5;
%{
% % ori(4049:4051,1)=0.5;
% % ori(4549:4551,1)=0.5;
% % ori(4749:4751,1)=0.5;
% % ori(4750,1)=1;
% % ori(5049:5051,1)=0.5;
% % ori(8229:8231,1)=0.5;
% % ori(8230,1)=1;
% % ori(4689:4691,1)=0.5;
% % ori(4690,1)=1;
% % ori(2046:2048,1)=0.5;
% % ori(2047,1)=1;
% % ori(2122:2124,1)=0.5;
% % ori(2123,1)=1;
% % ori(8330:8332,1)=0.5;
% % ori(8331,1)=1;
% % ori(4010:4012,1)=0.5;
% % ori(4011,1)=1;
%}
% lots of ground truth points can choose
pot_source = K * ori; %electrode potentials
pot_noise = en*randn(size(K,1),1); %random noise?

rms_signal = sqrt(sum(pot_source.^2)/length(pot_source));
rms_noise = sqrt(sum(pot_noise.^2)/length(pot_noise));

k = 10^(snr/20)*rms_noise/rms_signal; %the ratio
ori = k*ori; %degraded by the SNR
% Calculate the forward problem
pot = K * ori + pot_noise;
%% compute solution
fprintf('L1:    ');
pm.tol_rel = 1*10^(-3);
pm.MAX_ITERS = 1500;
pm.alpha2 = 100;%%%higher more sparse?
pm.rho = 0.2*pm.alpha2;%the ratio between rho and alpha is the smaller the quicker
[u, v, vt, history] = yl_l1(LFM_nrm, pot, pm, ss, 1); 
%% 3D display
figure(1);% show the result
title('Result');
ft_plot_vol(vol, 'edgecolor', 'none', 'facealpha', 0);
hold on;
Grid=sourcemodelwei;
gridin=Grid.pos(Grid.inside,:);
umax=max(abs(u));
u1 = u > max(u(:))*.95;% threshold to show the center of each reconstructed sources
scatter3(gridin(:,1),gridin(:,2),gridin(:,3),'filled','cdata',abs(u(Grid.inside)));
colorbar;
colormap('jet')

figure(2);% show the ground truth
title('Ground Truth');
ft_plot_vol(vol, 'edgecolor', 'none', 'facealpha', 0);
hold on;
Grid=sourcemodelwei;
gridin=Grid.pos(Grid.inside,:);
scatter3(gridin(:,1),gridin(:,2),gridin(:,3),'filled','cdata',abs(ori(Grid.inside)));
colorbar;
colormap('jet')