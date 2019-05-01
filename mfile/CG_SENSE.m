function [Image,para] = CG_SENSE(Data,para)
%[Image,para] = CG_SENSE(Data,para)
disp('Showing progress...')

ifplot = para.setting.plot;
ifGPU  = para.setting.ifGPU;
b = Data.first_est;

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
    Data.kSpace        = gpuArray(Data.kSpace);
    b          = gpuArray(b);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);

    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).kb_density_comp = gpuArray(Data.N(i).kb_density_comp);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
    
end

fidelity = @(im) compute_fidelity_yt_new(im,Data);
a = b;
b = zeros(size(b),'like',b);
p = a;
r = a;

aha = sqrt(sum(abs(vec(conj(a).*a)).^2));
rhr_1 = aha;

for iter_no = 1:para.Recon.noi

    if mod(iter_no,10) == 1
        t1 = tic;
    end
    rhr = rhr_1;
    delta(iter_no) = rhr/aha;

    q = -fidelity(p);
    phq = sqrt(sum(abs(vec(conj(p).*q)).^2));
    
    b(:,:,iter_no+1) = b(:,:,iter_no) + rhr/phq*p;
    r(:,:,iter_no+1) = r(:,:,iter_no) - rhr/phq*q;
    rhr_1 = sqrt(sum(abs(vec(conj(r(:,:,iter_no+1)).*r(:,:,iter_no+1))).^2));
    p = r(:,:,iter_no+1) + rhr_1/rhr*p;

    if ifplot == 1
        showImage(b(:,:,end),delta)
    end
    
    if mod(iter_no,10) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        toc(t1);
    end

end
para.delta = delta;
Image = flipud(gather(b(:,:,2:end)));

sx = size(Image,1);
filter = false(sx);
filter(sx/2+1,sx/2+1) = true;
filter = bwdist(filter);
filter = atan(100*(sx/2-filter)/(sx/2))/pi+0.5;
filter = fftshift2(filter);
Image = ifft2(filter.*fft2(Image));