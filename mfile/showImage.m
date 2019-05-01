function showImage(Image,Cost)
Image = (abs(Image));


figure(100);
clf;
subplot(2,2,[1 3])
imagesc(Image)
colormap gray
brighten(0.4)
axis image
axis off
        
subplot(2,2,[2 4])

hold on;
plot(log10(Cost));
legend('Total Cost')
xlabel 'iteration number'
ylabel 'log_{10}(delta)'

drawnow
end
