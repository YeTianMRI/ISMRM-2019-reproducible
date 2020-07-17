%%% For ISMRM 2019 reproducible research study group challenge
%%% University of Utah
%%% Utah Center for Advanced Imaging Research (UCAIR)
%%% The cardiovascular MRI group: Ganesh Adluru, Edward DiBella,
%%% Jason Mendes, Ye Tian and Brent Wilson
%%% contact: Ye Tian, ye.tian@utah.edu

% The code provided here is part of our radial reconstruction framework.
% The NUFFT codes were originally adopted from Jeffrey A. Fessler and
% were modulated by the authors to be GPU compatible. The ESPIRiT
% sensitivity map estimation requires the ESPIRiT package provided by
% Michael Lustig. To set up, run the 'setPath.m' function in the ESPIRiT
% package: run('ESPIRiT/setPath.m') and change the setting 
% 'para.setting.sens_map' be 'ESPIRIT'. This demo can also be run on GPU 
% by changing the 'para.setting.ifGPU' setting be '1'. The reconstruction 
% process can be visualized by setting 'para.setting.plot' be '1'.

% Our radial reconstruction frameworks can be found at:
% https://github.com/edibella/Reconstruction
% https://github.com/gadluru/Multiview-myocardial-perfusion-with-radial-SMS

% The original NUFFT code can be found at:
% https://web.eecs.umich.edu/~fessler/

% The ESPIRiT package can be downloaded at:
% https://people.eecs.berkeley.edu/~mlustig/Software.html

%% set up
clear
clc
close all
addpath mfile

para.setting.plot = 0;
para.setting.ifGPU = 0;
para.setting.sens_map = 'adaptive';

%% Load Brain data
rawdata_real = h5read('rawdata_brain_radial_96proj_12ch.h5','/rawdata');
trajectory   = h5read('rawdata_brain_radial_96proj_12ch.h5','/trajectory');
k_size = size(trajectory, 2);
im_size = max(trajectory(:)) * 2;

% modify data format to fit our reconstruction pipeline
kSpace = rawdata_real.r+1i*rawdata_real.i;
clear rawdata_real;
kSpace = permute(kSpace,[3,2,4,1]);
kSpace = kSpace*10^8; % scale up
trajectory = permute(trajectory,[2,1,3]);
kx = trajectory(:,:,1); kx = kx / im_size * k_size;
ky = trajectory(:,:,2); ky = ky / im_size * k_size;
clear trajectory

%% Reconstruction of reference
% initialize NUFFT
Data.kSpace = kSpace;
Data.N = NUFFT.init(kx,ky,2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);
Data.sens_map = get_sens_map(Data.first_est,para.setting.sens_map);

% coil combination
Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

% parameters
para.Recon.noi = 200;
Image_96_rays = CG_SENSE(Data,para);
Image_96_rays = crop_half_FOV(Image_96_rays, [im_size, im_size]);
Image_reference = Image_96_rays(:,:,end);

%% Conjugate Gradient NUFFT 48 rays
% initialize NUFFT
Data.kSpace = kSpace(:,1:2:end,:,:);
Data.N = NUFFT.init(kx(:,1:2:end),ky(:,1:2:end),2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);

% coil combination
Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

% parameters
para.Recon.noi = 100;
[Image_48_rays,para_48_rays] = CG_SENSE(Data,para);
Image_48_rays = crop_half_FOV(Image_48_rays, [im_size, im_size]);

%% Conjugate Gradient NUFFT 32 rays
% initialize NUFFT
Data.kSpace = kSpace(:,1:3:end,:,:);
Data.N = NUFFT.init(kx(:,1:3:end),ky(:,1:3:end),2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);

% coil combination
Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

% parameters
para.Recon.noi = 100;
[Image_32_rays,para_32_rays] = CG_SENSE(Data,para);
Image_32_rays = crop_half_FOV(Image_32_rays, [im_size, im_size]);

%% Conjugate Gradient NUFFT 24 rays
% initialize NUFFT
Data.kSpace = kSpace(:,1:4:end,:,:);
Data.N = NUFFT.init(kx(:,1:4:end),ky(:,1:4:end),2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);

% coil combination
Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

% parameters
para.Recon.noi = 100;
[Image_24_rays,para_24_rays] = CG_SENSE(Data,para);
Image_24_rays = crop_half_FOV(Image_24_rays, [im_size, im_size]);

%% Figure 4
difference_48_rays = abs(Image_48_rays - Image_reference)./abs(Image_reference);
difference_32_rays = abs(Image_32_rays - Image_reference)./abs(Image_reference);
difference_24_rays = abs(Image_24_rays - Image_reference)./abs(Image_reference);

difference_48_rays = log10(squeeze(sum(sum(difference_48_rays))));
difference_32_rays = log10(squeeze(sum(sum(difference_32_rays))));
difference_24_rays = log10(squeeze(sum(sum(difference_24_rays))));

figure
plot([difference_48_rays,difference_32_rays,difference_24_rays], 'LineWidth', 2)
xlabel 'number of iterations'
ylabel 'log_{10}(difference)'
legend('R = 2','R = 3','R = 4')
title 'Corresponding to Figure 4 of Ref. [1]'
set(gca, 'FontSize', 18)

%% single coil recon
Data.kSpace = kSpace(:,:,:,1);
Data.N = NUFFT.init(kx,ky,2,[6,6]);
Image_single_coil_96_rays = flipud(NUFFT.NUFFT_adj(Data.kSpace,Data.N));
Image_single_coil_96_rays = crop_half_FOV(Image_single_coil_96_rays, [im_size, im_size]);

Data.kSpace = kSpace(:,1:2:end,:,1);
Data.N = NUFFT.init(kx(:,1:2:end),ky(:,1:2:end),2,[6,6]);
Image_single_coil_48_rays = flipud(NUFFT.NUFFT_adj(Data.kSpace,Data.N));
Image_single_coil_48_rays = crop_half_FOV(Image_single_coil_48_rays, [im_size, im_size]);

Data.kSpace = kSpace(:,1:3:end,:,1);
Data.N = NUFFT.init(kx(:,1:3:end),ky(:,1:3:end),2,[6,6]);
Image_single_coil_32_rays = flipud(NUFFT.NUFFT_adj(Data.kSpace,Data.N));
Image_single_coil_32_rays = crop_half_FOV(Image_single_coil_32_rays, [im_size, im_size]);

Data.kSpace = kSpace(:,1:4:end,:,1);
Data.N = NUFFT.init(kx(:,1:4:end),ky(:,1:4:end),2,[6,6]);
Image_single_coil_24_rays = flipud(NUFFT.NUFFT_adj(Data.kSpace,Data.N));
Image_single_coil_24_rays = crop_half_FOV(Image_single_coil_24_rays, [im_size, im_size]);


%% Figure 5
figure
colormap gray
subplot(4,3,1)
imagesc(abs(Image_single_coil_96_rays)),axis image,axis off, title 'single coil R=1'
subplot(4,3,2)
imagesc(abs(Image_96_rays(:,:,1))),axis image,axis off, title 'Initial iteration = 1'
subplot(4,3,3)
imagesc(abs(Image_96_rays(:,:,3))),axis image,axis off, title 'Final iteration = 3'
subplot(4,3,4)
imagesc(abs(Image_single_coil_48_rays)),axis image,axis off, title 'single coil R=2'
subplot(4,3,5)
imagesc(abs(Image_48_rays(:,:,1))),axis image,axis off, title 'Initial iteration = 1'
subplot(4,3,6)
imagesc(abs(Image_48_rays(:,:,6))),axis image,axis off, title 'Final iteration = 6'
subplot(4,3,7)
imagesc(abs(Image_single_coil_32_rays)),axis image,axis off, title 'single coil R=3'
subplot(4,3,8)
imagesc(abs(Image_32_rays(:,:,1))),axis image,axis off, title 'Initial iteration = 1'
subplot(4,3,9)
imagesc(abs(Image_32_rays(:,:,15))),axis image,axis off, title 'Final iteration = 15'
subplot(4,3,10)
imagesc(abs(Image_single_coil_24_rays)),axis image,axis off, title 'single coil R=4'
subplot(4,3,11)
imagesc(abs(Image_24_rays(:,:,1))),axis image,axis off, title 'Initial iteration = 1'
subplot(4,3,12)
imagesc(abs(Image_24_rays(:,:,50))),axis image,axis off, title 'Final iteration = 50'

suptitle 'Corresponding to Figure 5 of Ref. [1]'
set(gca, 'FontSize', 18)

%% Load Cardiac Data
rawdata_real = h5read('rawdata_heart_radial_55proj_34ch.h5','/rawdata');
trajectory   = h5read('rawdata_heart_radial_55proj_34ch.h5','/trajectory');
k_size = size(trajectory, 2);
im_size = round(max(trajectory(:)) * 2);

% modify data format to fit our reconstruction pipeline
kSpace = rawdata_real.r+1i*rawdata_real.i;
clear rawdata_real;
kSpace = permute(kSpace,[3,2,4,1]);
kSpace = kSpace*10^8; % scale up
trajectory = permute(trajectory,[2,1,3]);
kx = trajectory(:,:,1); kx = kx / im_size * k_size;
ky = trajectory(:,:,2); ky = ky / im_size * k_size;
clear trajectory

para.Recon.noi = 100;
%% 55 projection
Data.kSpace = kSpace;
Data.N = NUFFT.init(kx,ky,2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);
Data.sens_map = get_sens_map(Data.first_est,para.setting.sens_map);

Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

Image_55_rays = CG_SENSE(Data,para);
Image_55_rays = abs(Image_55_rays(:,:,end));
Image_55_rays = crop_half_FOV(Image_55_rays, [im_size, im_size]);

%% 33 projection
Data.kSpace = kSpace(:,1:33,:,:);
Data.N = NUFFT.init(kx(:,1:33),ky(:,1:33),2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);

Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

Image_33_rays = CG_SENSE(Data,para);
Image_33_rays = abs(Image_33_rays(:,:,end));
Image_33_rays = crop_half_FOV(Image_33_rays, [im_size, im_size]);

%% 22 projection
Data.kSpace = kSpace(:,1:22,:,:);
Data.N = NUFFT.init(kx(:,1:22),ky(:,1:22),2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);

Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

Image_22_rays = CG_SENSE(Data,para);
Image_22_rays = abs(Image_22_rays(:,:,end));
Image_22_rays = crop_half_FOV(Image_22_rays, [im_size, im_size]);

%% 11 projection
Data.kSpace = kSpace(:,1:11,:,:);
Data.N = NUFFT.init(kx(:,1:11),ky(:,1:11),2,[6,6]);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace,Data.N);

Data.first_est = conj(Data.sens_map).*Data.first_est;
Data.first_est = sum(Data.first_est,4);

Image_11_rays = CG_SENSE(Data,para);
Image_11_rays = abs(Image_11_rays(:,:,end));
Image_11_rays = crop_half_FOV(Image_11_rays, [im_size, im_size]);

%% Figure 6
figure
colormap gray
subplot(1,4,1)
imagesc(Image_55_rays),axis image,axis off, title '55 projections'
subplot(1,4,2)
imagesc(Image_33_rays),axis image,axis off, title '33 projections'
subplot(1,4,3)
imagesc(Image_22_rays),axis image,axis off, title '22 projections'
subplot(1,4,4)
imagesc(Image_11_rays),axis image,axis off, title '11 projections'

suptitle 'Corresponding to Figure 6 of Ref. [1]'
set(gca, 'FontSize', 18)
