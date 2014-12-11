clear

%%Read in images, convert to double, grayscale and values between 0 to 1
pictures(:,:,1) = mat2gray(double(rgb2gray(imread('HalfCanon.jpg'))));
pictures(:,:,2) = mat2gray(double(rgb2gray(imread('HalfHolga.jpg'))));
pictures(:,:,3) = mat2gray(double(rgb2gray(imread('HalfScanner.jpg'))));
pictures(:,:,4) = mat2gray(double(rgb2gray(imread('HalfSony.jpg'))));

%%Construct a weight vector that looks like a V
SharpWeight = zeros (1, 512);
for i=1:255
    SharpWeight(256-i) = i/255;
    SharpWeight(257+i) = i/255;
end

%%Meshgrid and convert to polar coordinates.
N = 512;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

%Quantize to array with values between 1-100
QR = uint8(round(100*(R./max(R(:)))));
for i = 1:4
    %%Padarray, fft/fft2, and shift the images, also have values between 0 and 1
    ShiftedImage = fftshift(fft(padarray(sum(pictures(1:50, :, i)), [0 128])));
    ShiftedImage2 = fftshift(fft2(padarray(pictures(:,:,i), [128 128])));
    
    %Normalize images with the DC-component
    %PlotImage = ShiftedImage;
    %PlotImage2 = ShiftedImage2;
    ShiftedImage = abs(ShiftedImage./ShiftedImage(256));
    ShiftedImage2 = abs(ShiftedImage2./ShiftedImage2(256));
    
    %%Masking for second task in part one.
    for m = 1:100
        Maskm = (QR == m);
        masksum = sum(sum(Maskm));
        MaskedImage = ShiftedImage2(Maskm == 1);
        MaskedImageAvg(m) = mean(MaskedImage);
    end
    
    %Weight images and calculate sharpness
    Sharpness1(i) = sum(ShiftedImage.*SharpWeight);
    Sharpness2(i) = sum(MaskedImageAvg/max(MaskedImageAvg(:)).*(1:100));
end

%Plots for the labassignment for the Sony.
%{
subplot(3,1,1);
plot(real(PlotImage));
legend('Sony real part');
subplot(3,1,2);
plot(imag(PlotImage));
legend('Sony imaginary part');
subplot(3,1,3);
plot(abs(PlotImage));
legend('Sony absolute');
%}