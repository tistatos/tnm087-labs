clear all
load('winsuint8.mat')

%Create a meshgrid that is 64x64
N = 64;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

%Mask images with the grid, perform same operations as in part 1 on
%pictures, then take sum(sum(radial values .* pictures))
for m = 1:192
    pictures = fftshift(fft2(padarray(winsuint8(:,:,m), [16 16])));
    AbsPictures = abs(pictures/abs(pictures(33,33)));
    Sharpness(m) = sum(sum(R.*AbsPictures));
end

%Plot values
plot(Sharpness)
