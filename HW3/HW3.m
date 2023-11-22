%% HW3
clear; close all; clc; warning off;

%% Problem 1 Partial Fourier
load Data_Assignment3_Problem1.mat

% a) Partial fourier with factor of 5/8
img_dim = height(kspaceData_SingleCoil); 
kspace_partial_fourier = zeros(size(kspaceData_SingleCoil));
kspace_partial_fourier(1:img_dim*5/8,:) = kspaceData_SingleCoil(1:img_dim*5/8,:);
diff_img = zeros(size(kspaceData_SingleCoil));
diff_img(img_dim*5/8+1:end,:) = kspaceData_SingleCoil(img_dim*5/8+1:end,:);
figure 
subplot(3,2,1); imagesc(abs(ifftshift(ifft2(kspaceData_SingleCoil))))
title('Fully-sampled image (magnitude)'); axis equal; axis off; 
subplot(3,2,2); imagesc(angle(ifftshift(ifft2(kspaceData_SingleCoil))))
title('Fully-sampled image (phase)'); axis equal; axis off; 
subplot(3,2,3); imagesc(abs(ifftshift(ifft2(kspace_partial_fourier))))
title('Partial fourier image (magnitude)'); axis equal; axis off; 
subplot(3,2,4); imagesc(angle(ifftshift(ifft2(kspace_partial_fourier))))
title('Partial fourier image (phase)'); axis equal; axis off; 
subplot(3,2,5); imagesc(abs(ifftshift(ifft2(diff_img))))
title('Difference image (magnitude)'); axis equal; axis off; 
subplot(3,2,6); imagesc(angle(ifftshift(ifft2(diff_img))))
title('Difference image (phase)'); axis equal; axis off; colormap gray

% b) POCS recon
num_iterations = 50; 
phase_est_line_ind = [img_dim-img_dim*5/8+1:img_dim*5/8];
center_kspace = zeros(size(kspace_partial_fourier));
center_kspace(phase_est_line_ind,:) = kspace_partial_fourier(phase_est_line_ind,:);
center_phase = angle(fftshift(ifft2(center_kspace))); % phase image used to replace in the loop
% figure; imagesc(center_phase)
kspace_update = kspace_partial_fourier; % first kspace
for iter = 1:num_iterations
    % Get mag and phase from kspace_update (img domain)
    mag = abs(ifftshift(ifft2(kspace_update)));
    phase = angle(ifftshift(ifft2(kspace_update)));
%     figure; subplot(1,2,1); imagesc(mag); subplot(1,2,2); imagesc(phase)
    previous_kspace = kspace_update;
    
    % Replace phase with center_phase (img domain)
    phase = center_phase;
    
    % transform back to kspace (fourier domain)
    kspace_update = fft2(fftshift(mag.*exp(1i*phase)));
%     figure; imagesc(abs(kspace_update)); 
    
    % Replace this new kspace with original undersampled kspace (fourier domain)
    kspace_update(1:img_dim*5/8,:) = kspace_partial_fourier(1:img_dim*5/8,:);
%     figure; imagesc(abs(kspace_update)); 
    
end

figure
subplot(3,2,1); imagesc(abs(fftshift(ifft2(kspace_partial_fourier)))); axis image; axis off
title('Partial fourier image (mag)')
subplot(3,2,2); imagesc(angle(fftshift(ifft2(kspace_partial_fourier)))); axis image; axis off
title('Partial fourier image (phase)')
subplot(3,2,3); imagesc(abs(fftshift(ifft2(kspace_update)))); axis image; axis off
title('POCS recon after 50 iterations (mag)')
subplot(3,2,4); imagesc(angle(fftshift(ifft2(kspace_update)))); axis image; axis off
title('POCS recon after 50 iterations (phase)')
subplot(3,2,5); imagesc(abs(fftshift(ifft2(kspace_update-kspace_partial_fourier)))); axis image; axis off
title('Difference after 50 iterations (mag)')
subplot(3,2,6); imagesc(angle(fftshift(ifft2(kspace_update-kspace_partial_fourier)))); axis image; axis off; colormap gray
title('Difference after 50 iterations (phase)')

%% Problem 2 SENSE
load Data_Assignment3_Problem2.mat

% a) Coil-combined image
combine_img = zeros(size(kspaceData,1),size(kspaceData,1));
num_coil = size(coilmaps,3);
for coil = 1:num_coil
    combine_img = combine_img + conj(coilmaps(:,:,coil)).*ifftshift(ifft2(kspaceData(:,:,coil)));
end
figure; imagesc(abs(combine_img)); title('Coil-combined image (magnitude)'); axis image; axis off; colormap gray

% b) Aliased R=2 image
R2_kspace = zeros(size(kspaceData));
R2_kspace([1:2:end],:,:) = kspaceData([1:2:end],:,:);
% figure; subplot(1,2,1); imagesc(abs(kspaceData(:,:,3))); subplot(1,2,2); imagesc(abs(R2_kspace(:,:,3)))
combine_img_R2 = zeros(size(kspaceData,1),size(kspaceData,1));
for coil = 1:num_coil
    combine_img_R2 = combine_img_R2 + conj(coilmaps(:,:,coil)).*ifftshift(ifft2(R2_kspace(:,:,coil)));
end
combine_img_R2 = combine_img_R2./sqrt(sum(abs(coilmaps).^2,3)); % normalize
figure; imagesc(abs(combine_img_R2)); title('Image with acceleration R=2 (magnitude)'); axis image; axis off; colormap gray

% c) SENSE R=2 Reconstruction
R = 2; % acceleration factor
rho = SENSE(R2_kspace,coilmaps,R);
figure; subplot(1,2,1); imagesc(abs(rho)); axis image; axis off
title('SENSE recon R=2')
subplot(1,2,2); imagesc(abs(combine_img-rho)); axis image; axis off; colormap gray
title('Difference between fully sampled and SENSE recon R=2')

% c) SENSE R=4 Reconstruction
R = 4; % acceleration factor
R4_kspace = zeros(size(kspaceData));
R4_kspace([1:R:end],:,:) = kspaceData([1:R:end],:,:);
rho = SENSE(R4_kspace,coilmaps,R);
figure; subplot(1,2,1); imagesc(abs(rho)); axis image; axis off
title('SENSE recon R=4')
subplot(1,2,2); imagesc(abs(combine_img-rho)); axis image; axis off; colormap gray
title('Difference between fully sampled and SENSE recon R=4')

%% Function
function [rho] = SENSE(kspaceData,coilmaps,acceleration)
R = acceleration; % acceleration factor
% kspace = kspaceData([1:R:end],:,:); % reduced kspace matrix based on R
kspace = kspaceData;
mat_dim_long = length(kspaceData); % should be 200
mat_dim_short = size(kspace,1)/R; % should be 200/R 
num_coil = size(coilmaps,3);

% Get image in time domain
img = zeros(mat_dim_long,mat_dim_long,num_coil);
for coil = 1:num_coil
    img(:,:,coil) = ifftshift(ifft2(kspace(:,:,coil)));
end

% Build system of linear equation
rho = zeros(size(kspaceData,1),size(kspaceData,2)); % preallocate for final image
for col = 1:mat_dim_long
    for row = 1:mat_dim_short
        % for each pixel, get I from each coilmap. 
        % I should have length of 8 because there are 8 coils.
        I = squeeze(img(row,col,:)); 
        % C should be 8xR matrix
        C = transpose(squeeze(coilmaps([row:mat_dim_long/R:mat_dim_long],col,:)));
        % Solve for rho and populate image
        rho([row:mat_dim_long/R:mat_dim_long],col) = (C\I)*R;
    end
end
end