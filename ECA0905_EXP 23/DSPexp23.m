clc;
clear;
close all;

%% DSP IMAGE ENHANCEMENT, COMPRESSION, AND EDGE DETECTION

%% Input parameters
imgName = 'cameraman.tif';
medSize = [3 3];
jpegQuality = 50;
cannyTh = [0.1 0.2];
MAXI = 255;                         % Peak value for an 8-bit image

%% 1. Load image
I = imread(imgName);

if size(I, 3) == 3
    I = rgb2gray(I);
end

[rows, cols] = size(I);

fprintf('Image: %s | Size: %d x %d | Class: %s\n', ...
    imgName, rows, cols, class(I));
fprintf('Min = %d, Max = %d, Mean = %.2f, Std = %.2f\n', ...
    min(I(:)), max(I(:)), mean(double(I(:))), std(double(I(:))));

%% 2. Spatial filtering and enhancement
I_med = medfilt2(I, medSize);
I_enh = imadjust(I_med, stretchlim(I_med), []);

%% 3. JPEG compression and reconstruction
jpgFile = fullfile(tempdir, 'cameraman_q50.jpg');
imwrite(I, jpgFile, 'jpg', 'Quality', jpegQuality);
I_rec = imread(jpgFile);

origBytes = rows * cols;            % 8-bit grayscale: 1 byte per pixel
fileInfo = dir(jpgFile);
compBytes = fileInfo.bytes;

CR = origBytes / compBytes;
bpp = 8 * compBytes / (rows * cols);

%% 4. Edge detection
E_orig = edge(I, 'canny', cannyTh);
E_enh = edge(I_enh, 'canny', cannyTh);
E_rec = edge(I_rec, 'canny', cannyTh);

%% 5. Quality metrics
MSE = @(A, B) mean((double(A(:)) - double(B(:))).^2);

PSNR = @(A, B) 10 * log10(MAXI^2 / max(MSE(A, B), eps));

MSE_med = MSE(I, I_med);
PSNR_med = PSNR(I, I_med);

MSE_enh = MSE(I, I_enh);
PSNR_enh = PSNR(I, I_enh);

MSE_rec = MSE(I, I_rec);
PSNR_rec = PSNR(I, I_rec);

NIE_med = MSE_med / MAXI^2;
NIE_enh = MSE_enh / MAXI^2;
NIE_rec = MSE_rec / MAXI^2;

SSIM_rec = ssim(I_rec, I);

edgePix = [nnz(E_orig), nnz(E_enh), nnz(E_rec)];

%% 6. Display results
figure('Name', 'Enhancement, Compression, Edge Detection', ...
    'Color', 'w', 'Position', [40 40 1200 750]);

subplot(2,3,1);
imshow(I);
title('Original image');

subplot(2,3,2);
imshow(I_med);
title('After 3x3 median filter');

subplot(2,3,3);
imshow(I_enh);
title('Enhanced (median + imadjust)');

subplot(2,3,4);
imshow(I_rec);
title(sprintf('JPEG Q=%d reconstructed (CR = %.2f:1)', ...
    jpegQuality, CR));

subplot(2,3,5);
imshow(E_orig);
title('Canny edge map (original)');

subplot(2,3,6);
imshow(E_rec);
title('Canny edge map (JPEG reconstructed)');

% Histograms
figure('Name', 'Histograms', 'Color', 'w', ...
    'Position', [80 80 900 350]);

subplot(1,2,1);
imhist(I);
title('Histogram - original');

subplot(1,2,2);
imhist(I_enh);
title('Histogram - enhanced');

% Compression error and enhanced-image edges
figure('Name', 'Compression error', 'Color', 'w', ...
    'Position', [120 120 900 380]);

subplot(1,2,1);
imshow(abs(double(I) - double(I_rec)), []);
colorbar;
title('|Original - JPEG reconstructed|');

subplot(1,2,2);
imshow(E_enh);
title('Canny edge map (enhanced image)');

%% JPEG rate-distortion comparison
qList = 10:10:90;
psnrQ = zeros(size(qList));
crQ = zeros(size(qList));

for k = 1:numel(qList)
    tempFile = fullfile(tempdir, 'tmp_q.jpg');
    imwrite(I, tempFile, 'jpg', 'Quality', qList(k));

    reconstructed = imread(tempFile);
    tempInfo = dir(tempFile);

    psnrQ(k) = PSNR(I, reconstructed);
    crQ(k) = origBytes / tempInfo.bytes;
end

figure('Name', 'Rate-distortion', 'Color', 'w', ...
    'Position', [160 160 900 380]);

subplot(1,2,1);
plot(qList, psnrQ, '-o', 'LineWidth', 1.4);
grid on;
xlabel('JPEG quality');
ylabel('PSNR (dB)');
title('PSNR vs JPEG quality');

subplot(1,2,2);
plot(qList, crQ, '-s', 'LineWidth', 1.4);
grid on;
xlabel('JPEG quality');
ylabel('Compression ratio');
title('Compression ratio vs JPEG quality');

%% 7. Output parameter table
fprintf('\n=============== OUTPUT PARAMETER TABLE ===============\n');
fprintf('%-28s %10s %10s %12s\n', ...
    'Image', 'MSE', 'PSNR (dB)', 'MSE/MAXI^2');

fprintf('%-28s %10.2f %10.2f %12.6f\n', ...
    'Median filtered', MSE_med, PSNR_med, NIE_med);

fprintf('%-28s %10.2f %10.2f %12.6f\n', ...
    'Enhanced (imadjust)', MSE_enh, PSNR_enh, NIE_enh);

fprintf('%-28s %10.2f %10.2f %12.6f\n', ...
    'JPEG Q=50 reconstructed', MSE_rec, PSNR_rec, NIE_rec);

fprintf('\n--- Compression results (JPEG Q = %d) ---\n', jpegQuality);
fprintf('Original size  : %d bytes\n', origBytes);
fprintf('Compressed size: %d bytes\n', compBytes);
fprintf('Compression ratio: %.2f:1\n', CR);
fprintf('Bits per pixel : %.3f bpp\n', bpp);
fprintf('SSIM (JPEG vs original): %.4f\n', SSIM_rec);

fprintf('\n--- Edge detection (Canny [%.1f %.1f]) ---\n', cannyTh(1), cannyTh(2));
fprintf('Edge pixels - original image       : %d\n', edgePix(1));
fprintf('Edge pixels - enhanced image       : %d\n', edgePix(2));
fprintf('Edge pixels - JPEG reconstructed   : %d\n', edgePix(3));