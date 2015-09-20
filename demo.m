%This code demonstrate the reconstruction for GLIM using 3D deconvolution.
%Version 1: no sparse assumption added
clc;
clear all;
close all;
%Load the data
datafolder = 'E:\Data_for_embryo_tomography\Embryo_darta_0_55\';
measdatafile = strcat(datafolder,'Embryo_darta_0_55_zds_4x_16bit.tif');
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
imagesc(measdata(:,:,83));colormap gray;

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
psfdatapad=fftshift(circshift(psfdatapad,shiftamount));
psfdatapad = fftshift(fftn(psfdatapad));
measdata = fftshift(fftn(measdata));
figure(2);
subplot(121);
imagesc(log10(abs(measdata(:,:,end/2))));colorbar;
title('FFT of the data');
subplot(122);
imagesc(log10(abs(psfdatapad(:,:,end/2))));colorbar;
title('FFT of the derivative kernel')
outputdata = measdata.*conj(psfdatapad)./(abs(psfdatapad).^2+0.001);
clear measdata;
clear psfdatapad;
outputdata=real(ifftn(ifftshift(outputdata)));
figure(1);
subplot(133);
imagesc(outputdata(:,:,83));drawnow;colormap gray;

