function sens_map = get_sens_map(im,options)


switch options
    case 'adaptive'
        smooth = 10;
        im_for_sens = squeeze(sum(im,3));
        sens_map(:,:,1,:) = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens,smooth);
    case 'ESPIRIT'
        eigThresh_k = 0.02; % threshold of eigenvectors in k-space
        eigThresh_im = 0.9; % threshold of eigenvectors in image space
        kernel_size = [6,6];
        center_rad = 10;
        [sx,sy,~,coils,slice] = size(im);
        for i=1:slice
            im_for_sens = squeeze(mean(im(:,:,:,:,i),3));
            im_for_sens = fftshift(fftshift(im_for_sens,1),2);
            k_for_sens = fft2(im_for_sens);
            k_for_sens = fftshift(fftshift(k_for_sens,1),2);
            k_center = round(size(k_for_sens,1)/2);
            k_for_sens = k_for_sens(k_center-center_rad+1:k_center+center_rad,k_center-center_rad+1:k_center+center_rad,:);
            
            [k,S] = dat2Kernel(k_for_sens,kernel_size);
            idx = find(S >= S(1)*eigThresh_k,1,'last');
            [M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

            weights = (W - eigThresh_im)./(1-eigThresh_im).* (W> eigThresh_im);
            weights = -cos(pi*weights)/2 + 1/2;
            weights(weights==0) = 0.01;
            sens_map_temp = weights.*permute(M,[1,2,4,3]);
            sens_map_temp = sens_map_temp(:,:,end,:);
            sens_map(:,:,:,:,i) = permute(sens_map_temp,[1 2 5 4 3]);
        end

end

sens_map_scale = max(abs(sens_map(:)));
sens_map = sens_map/sens_map_scale;
sens_map_conj = conj(sens_map);

sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);

sens_correct_term = sqrt(sens_correct_term);
sens_map = bsxfun(@times,sens_correct_term,sens_map);

sens_map = single(sens_map);