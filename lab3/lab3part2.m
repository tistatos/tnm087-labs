clear all
load('winsuint8.mat')

N = 64;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

for m = 1:192
    
    pictures = fftshift(fft2(padarray(winsuint8(:,:,m), [16 16])));
    AbsPictures = abs(pictures/abs(pictures(33,33)));
 
    SR = R ./ R(N/2 - 1, 1);
    Sharpness(m) = sum(sum(SR.*AbsPictures));
    
end


plot(Sharpness)
