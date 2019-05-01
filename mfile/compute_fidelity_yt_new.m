function [fidelity_update,fidelity_norm] = compute_fidelity_yt_new(image,Data)

fidelity_update = bsxfun(@times,image,Data.sens_map);
fidelity_update = NUFFT.NUFFT(fidelity_update,Data.N);

fidelity_update = Data.kSpace - fidelity_update;
fidelity_norm = sum(abs(fidelity_update(:)).^2);

fidelity_update = NUFFT.NUFFT_adj(fidelity_update,Data.N);
fidelity_update = bsxfun(@times,fidelity_update,Data.sens_map_conj);
fidelity_update = sum(fidelity_update,4);

end
