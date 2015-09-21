%This code demonstrate the reconstruction for GLIM using 3D deconvolution.
%Version 1: no sparse assumption added
clc;
clear all;
close all;
%Load the data
%datafolder = 'E:\Data_for_embryo_tomography\Embryo_darta_0_55\';
datafolder = 'E:\Data_for_embryo_tomography\Four_half_um_beads\'
%measdatafile = strcat(datafolder,'Embryo_darta_0_55_zds_4x_16bit.tif');
measdatafile = strcat(datafolder,'beads.tif');
measdatainfo = imfinfo(measdatafile);
nz = length(measdatainfo);
frame1 = imread(measdatafile,1);
nr = size(frame1,1);
nc = size(frame1,2);
measdata = zeros(nr,nc,nz);
for zidx = 1:nz
    disp(['Reading frame: ' num2str(zidx)]);
    measdata(:,:,zidx)=imread(measdatafile,zidx);    
end

ds_fact = 2;
measdata = measdata(1:ds_fact:end,1:ds_fact:end,:); %Downsample by a factor of 2 for more r
figure(1);
subplot(131);
imagesc(measdata(:,:,round(end/2)));colormap gray;

%Read the PSF
psfdatafile = strcat(datafolder,'PSF_data_0_55.tif');
psfdatainfo = imfinfo(psfdatafile);
nzpsf = length(psfdatainfo);
psfframe1 = imread(psfdatafile,1);
nrpsf = size(psfframe1,1);
ncpsf = size(psfframe1,2);
psfdata = zeros(nrpsf,ncpsf,nzpsf);
for zidxpsf = 1:nzpsf
    disp(['PSF reading frame: ', num2str(zidxpsf)]);
    psfdata(:,:,zidxpsf) = imread(psfdatafile,zidxpsf);
end
psfdata = psfdata(1:ds_fact:end,1:ds_fact:end,:);
psfdata = psfdata/sum(psfdata(:));
%measdata = imread
figure(1);
subplot(132);
imagesc(psfdata(:,:,5));colormap gray;
center_coord = [round(61/ds_fact) round(67/ds_fact) 3];%[y,x,z] %Coordinates of the central frequency to perform the shiftting
shiftamount = -center_coord;
psfdatapad = zeros(size(measdata));
psfdatapad(1:size(psfdata,1),1:size(psfdata,2),1:size(psfdata,3))=psfdata;
%Shift the data so that the psf is centered at the corner
psfdatapad=circshift(psfdatapad,shiftamount);
nxds = size(psfdatapad,2);%Number of samples in x after downsampling
nyds = size(psfdatapad,1);
nzds = size(psfdatapad,3);
[kx_arr,ky_arr,kz_arr]=meshgrid(linspace(-nxds/2,nxds/2,nxds),linspace(-nyds/2,nyds/2,nyds),linspace(-nzds/2,nzds/2,nzds));
%defining the bandwidth for (kx,ky,kz)
km =80; 
mask = ((kx_arr.^2 + ky_arr.^2 + kz_arr.^2)<km^2); %Maximum spatial frequency that can be covered. May need a tapper here..
clear kx_arr;
clear ky_arr;
clear kz_arr;
psfdatapad =ifftshift(fftshift(fftn(psfdatapad)).*mask);
measdata = (fftn(measdata));
figure(2);
subplot(121);
imagesc(log10(abs(measdata(:,:,round(end/2)))));colorbar;
title('FFT of the data');
subplot(122);
imagesc(log10(abs(psfdatapad(:,:,round(end/2)))));colorbar;
title('FFT of the derivative kernel')
outputdata = measdata.*conj(psfdatapad)./(abs(psfdatapad).^2+5);
%outputdata = measdata;
clear measdata;
clear psfdatapad;
outputdata=real(ifftn(outputdata));
figure(1);
subplot(133);
imagesc(outputdata(:,:,round(end/2)));drawnow;colormap gray;
for zidx = 1:nzds
    disp(['Saving z = ' num2str(zidx)]);
    writeTIFF(cast(real(outputdata(:,:,zidx)),'single'),strcat(datafolder,'recon_z_',num2str(zidx),'.tif'));
end