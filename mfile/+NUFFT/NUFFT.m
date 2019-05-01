function kSpace_radial = NUFFT(image,N)

%sx = N.siz(1);
sx_over = N.sx_over;

%nof = N.siz(3);
%nc = size(image,4);
%nSMS = size(image,5);
%throw_away_x = (sx_over-sx)/2;
%throw_away_y = (sx_over-sx)/2;
%throw_away_x = round(throw_away_x);
%throw_away_y = round(throw_away_y);
%image(1:2:end,:,:,:) = -image(1:2:end,:,:,:);
%image(:,1:2:end,:,:) = -image(:,1:2:end,:,:);

%image_over = zeros(sx_over,sx_over,nof,nc,nSMS);

image = bsxfun(@times,image,N.kb_density_comp);
%image_over(1:sx,1:sx,:,:,:) = image;

%image_over(1:2:end,:,:,:) = - image_over(1:2:end,:,:,:);
%image_over(:,1:2:end,:,:) = - image_over(:,1:2:end,:,:);

%image_over = fftshift(image_over,1);
%image_over = fftshift(image_over,2);

%for i=1:3
%kSpace_cart(:,:,:,:,i) = fft2(image(:,:,:,:,i),sx_over,sx_over);
%end
if N.sx_over < N.siz(1)
    image = fftshift2(image);
end
kSpace_cart = fft2(image,sx_over,sx_over);

%kSpace_cart = fftshift(kSpace_cart,1);
%kSpace_cart = fftshift(kSpace_cart,2);
%kSpace_cart(1:2:end,:,:,:) = - kSpace_cart(1:2:end,:,:,:);
%kSpace_cart(:,1:2:end,:,:) = - kSpace_cart(:,1:2:end,:,:);

kSpace_radial = NUFFT.cart2rad(kSpace_cart,N);