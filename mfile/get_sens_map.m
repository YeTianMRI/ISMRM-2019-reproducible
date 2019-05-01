function sens_map = get_sens_map(im)

smooth = 10;

im_for_sens = squeeze(sum(im,3));

sens_map(:,:,1,:) = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens,smooth);

sens_map_scale = max(abs(sens_map(:)));
sens_map = sens_map/sens_map_scale;
sens_map_conj = conj(sens_map);

sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);

sens_correct_term = sqrt(sens_correct_term);
sens_map = bsxfun(@times,sens_correct_term,sens_map);

sens_map = single(sens_map);