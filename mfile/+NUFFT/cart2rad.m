function kSpace_radial = cart2rad_new(kSpace_cart,N)

%sx = N.siz(1);
%nor = N.siz(2);
nof = N.siz(3);
nc = size(kSpace_cart,4);
nSMS = size(kSpace_cart,5);
sx_over = N.sx_over;

kSpace_cart = reshape(kSpace_cart,[sx_over*sx_over*nof,nc*nSMS]);

%kSpace_cart = permute(kSpace_cart,[1 3 2 4]);

%kSpace_radial = single(zeros([sx*nor*nof,nc,nSMS]));
%if isa(kSpace_cart,'gpuArray')
%    kSpace_radial = gpuArray(kSpace_radial);
%end
%for j=1:nSMS
%for i=1:nof
kSpace_radial = single(N.S'*double(kSpace_cart));
%end
%end

%kSpace_radial = reshape(kSpace_radial,[sx*nor,prod(core_size),nof,nc,nSMS]);
%kSpace_radial = permute(kSpace_radial,[1 3 2 4 5]);

%kSpace_radial = bsxfun(@times,kSpace_radial,N.weight_kb);

%kSpace_radial = sum(kSpace_radial,3);

kSpace_radial = reshape(kSpace_radial,[N.siz,nc,nSMS]);