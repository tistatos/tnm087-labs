% Lab 3 - Operations in the fourier domain
clear 


pictures(:,:,1) = double(rgb2gray(imread('HalfCanon.jpg')));
pictures(:,:,2) = double(rgb2gray(imread('HalfHolga.jpg')));
pictures(:,:,3) = double(rgb2gray(imread('HalfScanner.jpg')));
pictures(:,:,4) = double(rgb2gray(imread('HalfSony.jpg')));

HalfCanon = pictures(:,:,1)/max(max(pictures(:,:,1)));
HalfHolga = pictures(:,:,2)/max(max(pictures(:,:,2)));
HalfScanner = pictures(:,:,3)/max(max(pictures(:,:,3)));
HalfSony = pictures(:,:,4)/max(max(pictures(:,:,4)));

EdgeCanon = HalfCanon(1:50, :);
EdgeHolga = HalfHolga(1:50, :);
EdgeScanner = HalfScanner(1:50, :);
EdgeSony = HalfSony(1:50, :);

padCanon = padarray(sum(EdgeCanon), 127);
padHolga = padarray(sum(EdgeHolga), 127);
padScanner = padarray(sum(EdgeScanner), 127);
padSony = padarray(sum(EdgeSony), 127);

FFT1EdgeCanon = fftshift(fft(padCanon));
FFT1EdgeHolga = fftshift(fft(padHolga));
FFT1EdgeScanner = fftshift(fft(padScanner));
FFT1EdgeSony = fftshift(fft(padSony));

%plot(real(FFT1EdgeCanon), imag(FFT1EdgeCanon));
%plot(real(FFT1EdgeHolga), imag(FFT1EdgeHolga));
%plot(real(FFT1EdgeScanner), imag(FFT1EdgeScanner));
%plot(real(FFT1EdgeSony), imag(FFT1EdgeSony));

AbsFFT1EdgeCanon = abs(FFT1EdgeCanon/(length(FFT1EdgeCanon)/2));
AbsFFT1EdgeHolga = abs(FFT1EdgeHolga/(length(FFT1EdgeHolga)/2));
AbsFFT1EdgeScanner = abs(FFT1EdgeScanner/(length(FFT1EdgeScanner)/2));
AbsFFT1EdgeSony = abs(FFT1EdgeSony/(length(FFT1EdgeSony)/2));

%Sharpness test using standard deviation. Might work?
sharpness1 = sharpness_measure(AbsFFT1EdgeCanon);
sharpness2 = sharpness_measure(AbsFFT1EdgeHolga);
sharpness3 = sharpness_measure(AbsFFT1EdgeScanner);
sharpness4 = sharpness_measure(AbsFFT1EdgeSony);

T1 = table(sharpness1, sharpness2, sharpness3, sharpness4);

%% Part b

%For the meshgrid
N = 50;

EdgeCanon2 = padarray(HalfCanon, [128 128]);
EdgeHolga2 = padarray(HalfHolga, [128 128]);
EdgeScanner2 = padarray(HalfScanner, [128 128]);
EdgeSony2 = padarray(HalfSony, [128 128]);

EdgeCanon2 = fftshift(fft2(EdgeCanon2));
EdgeHolga2 = fftshift(fft2(EdgeHolga2));
EdgeScanner2 = fftshift(fft2(EdgeScanner2));
EdgeSony2 = fftshift(fft2(EdgeSony2));

AbsEdgeCanon2 = abs(EdgeCanon2/(length(EdgeCanon2)/2));
AbsEdgeHolga2 = abs(EdgeHolga2/(length(EdgeHolga2)/2));
AbsEdgeScanner2 = abs(EdgeScanner2/(length(EdgeScanner2)/2));
AbsEdgeSony2 = abs(EdgeSony2/(length(EdgeSony2)/2));

%plot(AbsEdgeCanon2, 'red');
%hold on
%plot(AbsEdgeHolga2, 'blue');
%hold on
%plot(AbsEdgeScanner2, 'green');
%hold on
%plot(AbsEdgeSony2, 'black');

[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

SR = R ./ R(N/2 - 1, 1);
QR = uint8(round(100*(SR./max(SR(:)))));


canonMaskAvg = zeros(100, 1);
HolgaMaskAvg = zeros(100, 1);
for m = 1:100
    Maskm = QR == m;
    %Sum over the pixel values
    masksum = sum(sum(Maskm));
    %object point?
    Canonmask = AbsEdgeCanon2(Maskm == 1);
    HolgaMask = AbsEdgeHolga2(Maskm == 1);
    canonMaskAvg(m) = mean(Canonmask);
    HolgaMaskAvg(m) = mean(HolgaMask);
end

plot(canonMaskAvg/max(canonMaskAvg)); 
hold on
plot(HolgaMaskAvg/max(HolgaMaskAvg), 'red'); 

%Sharpness test using standard deviation. Might work?
sharpness5 = sharpness_measure(AbsEdgeCanon2);
sharpness6 = sharpness_measure(AbsEdgeHolga2);
sharpness7 = sharpness_measure(AbsEdgeScanner2);
sharpness8 = sharpness_measure(AbsEdgeSony2);


T2 = table(sharpness5, sharpness6, sharpness7, sharpness8);