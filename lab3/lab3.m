% Lab 3 - Operations in the fourier domain
clear 

pictures(:,:,1) = double(rgb2gray(imread('HalfCanon.jpg')));
pictures(:,:,2) = double(rgb2gray(imread('HalfHolga.jpg')));
pictures(:,:,3) = double(rgb2gray(imread('HalfScanner.jpg')));
pictures(:,:,4) = double(rgb2gray(imread('HalfSony.jpg')));

SharpWeight = zeros (1, 512);
for i=1:255
    SharpWeight(256-i) = i/25;
    SharpWeight(257+i) = i/25;
end

for i = 1:4
    pictures(:,:,i) = pictures(:,:,i)/max(max(pictures(:,:,i)));
    EdgePicture = pictures(1:50, :, i);
    ShiftedImages = fftshift(fft(padarray(sum(EdgePicture), [0 128])));
    ShiftedImages = ShiftedImages/ShiftedImages(256);
    Sharpness(i) = sum(abs(ShiftedImages).*SharpWeight);
end

HalfCanon = pictures(:,:,1)/max(max(pictures(:,:,1)));
HalfHolga = pictures(:,:,2)/max(max(pictures(:,:,2)));
HalfScanner = pictures(:,:,3)/max(max(pictures(:,:,3)));
HalfSony = pictures(:,:,4)/max(max(pictures(:,:,4)));

EdgeCanon = HalfCanon(1:50, :);
EdgeHolga = HalfHolga(1:50, :);
EdgeScanner = HalfScanner(1:50, :);
EdgeSony = HalfSony(1:50, :);

padCanon = padarray(sum(EdgeCanon), [0 128]);
padHolga = padarray(sum(EdgeHolga), [0 128]);
padScanner = padarray(sum(EdgeScanner), [0 128]);
padSony = padarray(sum(EdgeSony), [0 128]);

FFT1EdgeCanon = fftshift(fft(padCanon));
FFT1EdgeHolga = fftshift(fft(padHolga));
FFT1EdgeScanner = fftshift(fft(padScanner));
FFT1EdgeSony = fftshift(fft(padSony));

%plot(real(FFT1EdgeCanon), imag(FFT1EdgeCanon));
%plot(real(FFT1EdgeHolga), imag(FFT1EdgeHolga));
%plot(real(FFT1EdgeScanner), imag(FFT1EdgeScanner));
%plot(real(FFT1EdgeSony), imag(FFT1EdgeSony));

AbsFFT1EdgeCanon = abs(FFT1EdgeCanon/FFT1EdgeCanon(256));
AbsFFT1EdgeHolga = abs(FFT1EdgeHolga/FFT1EdgeHolga(256));
AbsFFT1EdgeScanner = abs(FFT1EdgeScanner/FFT1EdgeScanner(256));
AbsFFT1EdgeSony = abs(FFT1EdgeSony/FFT1EdgeSony(256));

CanonSharp = sum(AbsFFT1EdgeCanon.*SharpWeight);
HolgaSharp = sum(AbsFFT1EdgeHolga.*SharpWeight);
ScannerSharp = sum(AbsFFT1EdgeScanner.*SharpWeight);
SonySharp = sum(AbsFFT1EdgeSony.*SharpWeight);

T1 = table(CanonSharp, HolgaSharp, ScannerSharp, SonySharp);

%Sharpness

%% Part b
EdgeCanon2 = padarray(HalfCanon, [128 128]);
EdgeHolga2 = padarray(HalfHolga, [128 128]);
EdgeScanner2 = padarray(HalfScanner, [128 128]);
EdgeSony2 = padarray(HalfSony, [128 128]);

EdgeCanon2 = fftshift(fft2(EdgeCanon2));
EdgeHolga2 = fftshift(fft2(EdgeHolga2));
EdgeScanner2 = fftshift(fft2(EdgeScanner2));
EdgeSony2 = fftshift(fft2(EdgeSony2));

AbsEdgeCanon2 = abs(EdgeCanon2/EdgeCanon2(256));
AbsEdgeHolga2 = abs(EdgeHolga2/EdgeHolga2(256));
AbsEdgeScanner2 = abs(EdgeScanner2/EdgeScanner2(256));
AbsEdgeSony2 = abs(EdgeSony2/EdgeSony(256));

%plot(AbsEdgeCanon2, 'red');
%hold on
%plot(AbsEdgeHolga2, 'blue');
%hold on
%plot(AbsEdgeScanner2, 'green');
%hold on
%plot(AbsEdgeSony2, 'black');

SharpWeight = zeros (1, 100);
for i=1:50
    SharpWeight(51-i) = i/5;
    SharpWeight(50+i) = i/5;
end

%Calculate sharpness
N = 512;
[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

SR = R ./ R(N/2 - 1, 1);

%Quantize to array with values between 1-100
QR = uint8(round(100*(SR./max(SR(:)))));


for m = 1:100  
    Maskm = (QR == m);
    masksum = sum(sum(Maskm));
    CanonMask = AbsEdgeCanon2(Maskm == 1);
    HolgaMask = AbsEdgeHolga2(Maskm == 1);
    ScannerMask = AbsEdgeScanner2(Maskm == 1);
    SonyMask = AbsEdgeSony2(Maskm == 1);

   

    CanonMaskAvg(m) = mean(CanonMask); 
    HolgaMaskAvg(m) = mean(HolgaMask);
    ScannerMaskAvg(m) = mean(ScannerMask);
    SonyMaskAvg(m) = mean(SonyMask);
end


CanonSharp = sum(CanonMaskAvg/max(CanonMaskAvg(:)).*SharpWeight)
HolgaSharp = sum(HolgaMaskAvg/max(HolgaMaskAvg(:)).*SharpWeight)
ScannerSharp = sum(ScannerMaskAvg/max(ScannerMaskAvg(:)).*SharpWeight)
SonySharp = sum(SonyMaskAvg/max(SonyMaskAvg(:)).*SharpWeight)


T2 = table(CanonSharp, HolgaSharp, ScannerSharp, SonySharp);