clc;
clear;
close all;
rng(1);

%% INPUT PARAMETERS
P.N = 64;
P.pilots = [8 22 44 58];          % 1-based pilot indices
P.Lcp = 16;
P.Xp = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
P.dataIdx = setdiff(1:P.N, P.pilots);

SNR_dB = 20;
NsymMain = 3000;
NsymSweep = 4000;
snrSweep = 0:5:30;

% Four-tap exponential power-delay profile, normalized to unit power
pdp = exp(-(0:3)/1.5);
pdp = pdp / sum(pdp);
P.pdp = pdp(:);

% Channel frequency-correlation matrix
k = (0:P.N-1).';
dk = k - k.';
P.R = zeros(P.N);

for ell = 0:numel(pdp)-1
    P.R = P.R + pdp(ell+1) .* ...
        exp(-1j*2*pi*dk*ell/P.N);
end

%% MAIN RUN AT 20 dB
out = run_ofdm(SNR_dB, NsymMain, P);

%% BER AND MSE SWEEP
berLS = zeros(size(snrSweep));
berMMSE = zeros(size(snrSweep));
berPerfect = zeros(size(snrSweep));
mseLS = zeros(size(snrSweep));
mseMMSE = zeros(size(snrSweep));

for ii = 1:numel(snrSweep)
    result = run_ofdm(snrSweep(ii), NsymSweep, P);
    berLS(ii) = result.ber(1);
    berMMSE(ii) = result.ber(2);
    berPerfect(ii) = result.ber(3);
    mseLS(ii) = result.mse(1);
    mseMMSE(ii) = result.mse(2);
end

% Theoretical QPSK BER in Rayleigh fading with ideal channel knowledge
gamma_b = 10.^(snrSweep/10) / 2;
berTheory = 0.5 * (1 - sqrt(gamma_b ./ (1 + gamma_b)));

%% FIGURE 1: CHANNEL ESTIMATES
k = 1:P.N;

figure('Name', 'Channel estimates', 'Color', 'w');
subplot(1,2,1);
plot(k, abs(out.H1(:,1)), 'k', 'LineWidth', 1.8);
hold on;
plot(k, abs(out.H1(:,2)), 'r--', 'LineWidth', 1.3);
plot(k, abs(out.H1(:,3)), 'b-.', 'LineWidth', 1.3);
plot(P.pilots, abs(out.Hpil1), 'go', 'MarkerFaceColor', 'g');
grid on;
xlabel('Subcarrier index');
ylabel('|H(k)|');
title(sprintf('Channel magnitude (SNR = %d dB)', SNR_dB));
legend('Actual', 'LS (interpolated)', 'MMSE', 'LS at pilots', ...
    'Location', 'best');

subplot(1,2,2);
plot(k, angle(out.H1(:,1)), 'k', 'LineWidth', 1.8);
hold on;
plot(k, angle(out.H1(:,2)), 'r--', 'LineWidth', 1.3);
plot(k, angle(out.H1(:,3)), 'b-.', 'LineWidth', 1.3);
grid on;
xlabel('Subcarrier index');
ylabel('Phase (rad)');
title('Channel phase');
legend('Actual', 'LS', 'MMSE', 'Location', 'best');

%% FIGURE 2: ESTIMATION ERROR
figure('Name', 'Estimation error', 'Color', 'w');
subplot(1,2,1);
semilogy(k, max(out.mse_k(:,1), eps), 'r-o', 'MarkerSize', 4);
hold on;
semilogy(k, max(out.mse_k(:,2), eps), 'b-s', 'MarkerSize', 4);
grid on;
xlabel('Subcarrier index');
ylabel('MSE = E|H - H_{est}|^2');
title(sprintf('Per-subcarrier estimation error (%d dB)', SNR_dB));
legend('LS', 'MMSE', 'Location', 'best');

subplot(1,2,2);
semilogy(snrSweep, max(mseLS, eps), 'r-o', 'LineWidth', 1.4);
hold on;
semilogy(snrSweep, max(mseMMSE, eps), 'b-s', 'LineWidth', 1.4);
grid on;
xlabel('SNR (dB)');
ylabel('Channel estimation MSE');
title('Estimation MSE vs SNR');
legend('LS', 'MMSE', 'Location', 'best');

%% FIGURE 3: EQUALIZED CONSTELLATIONS
figure('Name', 'Constellations', 'Color', 'w');
titles = {'Perfect CSI', 'LS estimate', 'MMSE estimate'};
fields = {'perfect', 'ls', 'mmse'};
berIndices = [3 1 2];
ideal = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2);

for ii = 1:3
    subplot(1,3,ii);
    z = out.const.(fields{ii});
    plot(real(z), imag(z), '.', 'MarkerSize', 3);
    hold on;
    plot(real(ideal), imag(ideal), 'r+', ...
        'MarkerSize', 12, 'LineWidth', 2);
    grid on;
    axis equal;
    axis([-2 2 -2 2]);
    xlabel('In-phase');
    ylabel('Quadrature');
    title(sprintf('%s (BER = %.2e)', ...
        titles{ii}, out.ber(berIndices(ii))));
end

%% FIGURE 4: BER VS SNR
figure('Name', 'BER comparison', 'Color', 'w');
semilogy(snrSweep, max(berLS, 1e-6), 'r-o', 'LineWidth', 1.5);
hold on;
semilogy(snrSweep, max(berMMSE, 1e-6), 'b-s', 'LineWidth', 1.5);
semilogy(snrSweep, max(berPerfect, 1e-6), 'k-^', 'LineWidth', 1.5);
semilogy(snrSweep, berTheory, 'g--', 'LineWidth', 1.2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit error rate');
title('BER: LS vs MMSE channel estimation');
ylim([1e-6 1]);
legend('LS', 'MMSE', 'Perfect CSI', ...
    'Theory (Rayleigh, ideal CSI)', 'Location', 'southwest');

%% OUTPUT TABLES
fprintf('\n=========== RESULTS AT SNR = %d dB ===========\n', SNR_dB);
fprintf('Channel estimates at pilot subcarriers, first OFDM symbol:\n');
fprintf('%-8s %-22s %-22s %-22s\n', ...
    'Pilot', 'Actual H', 'LS estimate', 'MMSE estimate');

for ii = 1:numel(P.pilots)
    pilot = P.pilots(ii);
    fprintf('%-8d %-22s %-22s %-22s\n', pilot, ...
        cstr(out.H1(pilot,1)), ...
        cstr(out.H1(pilot,2)), ...
        cstr(out.H1(pilot,3)));
end

fprintf('\n%-30s %-14s %-14s\n', 'Metric', 'LS', 'MMSE');
fprintf('%-30s %-14.4e %-14.4e\n', ...
    'Channel estimation MSE', out.mse(1), out.mse(2));
fprintf('%-30s %-14.4e %-14.4e\n', ...
    'Bit error rate (BER)', out.ber(1), out.ber(2));
fprintf('Perfect-CSI BER: %.4e\n', out.ber(3));
fprintf('Bits simulated at %d dB: %d\n', ...
    SNR_dB, 2*numel(P.dataIdx)*NsymMain);

fprintf('\n%-10s %-12s %-12s %-12s %-12s %-12s\n', ...
    'SNR (dB)', 'BER LS', 'BER MMSE', 'BER perfect', 'MSE LS', 'MSE MMSE');

for ii = 1:numel(snrSweep)
    fprintf('%-10d %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e\n', ...
        snrSweep(ii), berLS(ii), berMMSE(ii), berPerfect(ii), ...
        mseLS(ii), mseMMSE(ii));
end

%% LOCAL FUNCTIONS

function s = cstr(z)
    s = sprintf('%.3f%+.3fj', real(z), imag(z));
end

function out = run_ofdm(snr_dB, Nsym, P)

    N = P.N;
    pilots = P.pilots(:);
    Np = numel(pilots);
    dataIdx = P.dataIdx(:);
    Nd = numel(dataIdx);
    channelLength = numel(P.pdp);

    noiseVar = 10^(-snr_dB/10);

    %% Transmitter: QPSK data and known pilots
    b1 = randi([0 1], Nd, Nsym);
    b2 = randi([0 1], Nd, Nsym);
    dataSym = ((1-2*b1) + 1j*(1-2*b2)) / sqrt(2);

    X = zeros(N, Nsym);
    X(dataIdx,:) = dataSym;
    X(pilots,:) = repmat(P.Xp(:), 1, Nsym);

    x = sqrt(N) * ifft(X, [], 1);
    xcp = [x(N-P.Lcp+1:N,:); x];

    %% Rayleigh fading channel and AWGN
    h = sqrt(P.pdp/2) .* ...
        (randn(channelLength,Nsym) + 1j*randn(channelLength,Nsym));

    y = zeros(N + P.Lcp, Nsym);
    for ss = 1:Nsym
        y(:,ss) = filter(h(:,ss), 1, xcp(:,ss));
    end

    y = y + sqrt(noiseVar/2) .* ...
        (randn(size(y)) + 1j*randn(size(y)));

    %% Receiver: remove cyclic prefix and FFT
    Y = fft(y(P.Lcp+1:end,:), [], 1) / sqrt(N);
    Hactual = fft(h, N, 1);

    %% LS channel estimation at pilot positions
    HlsPilot = Y(pilots,:) ./ P.Xp(:);
    Hls = interp1(pilots, HlsPilot, (1:N).', 'linear', 'extrap');

    %% MMSE channel estimation
    Rpp = P.R(pilots,pilots);
    Rfp = P.R(:,pilots);
    W = Rfp / (Rpp + noiseVar*eye(Np));
    Hmmse = W * HlsPilot;

    %% One-tap zero-forcing equalization
    Ydata = Y(dataIdx,:);
    Xperfect = Ydata ./ Hactual(dataIdx,:);
    Xls = Ydata ./ Hls(dataIdx,:);
    Xmmse = Ydata ./ Hmmse(dataIdx,:);

    %% BER in fixed order: LS, MMSE, perfect CSI
    berLS = qpsk_ber(Xls, b1, b2);
    berMMSE = qpsk_ber(Xmmse, b1, b2);
    berPerfect = qpsk_ber(Xperfect, b1, b2);
    out.ber = [berLS, berMMSE, berPerfect];

    %% Channel estimation MSE
    errLS = abs(Hactual - Hls).^2;
    errMMSE = abs(Hactual - Hmmse).^2;

    out.mse = [mean(errLS(:)), mean(errMMSE(:))];
    out.mse_k = [mean(errLS,2), mean(errMMSE,2)];

    %% Values for plots
    out.H1 = [Hactual(:,1), Hls(:,1), Hmmse(:,1)];
    out.Hpil1 = HlsPilot(:,1);

    nPlot = min(300, Nsym);
    out.const.perfect = reshape(Xperfect(:,1:nPlot), [], 1);
    out.const.ls = reshape(Xls(:,1:nPlot), [], 1);
    out.const.mmse = reshape(Xmmse(:,1:nPlot), [], 1);
end

function ber = qpsk_ber(z, b1, b2)
    errI = (real(z) < 0) ~= (b1 == 1);
    errQ = (imag(z) < 0) ~= (b2 == 1);
    ber = (sum(errI(:)) + sum(errQ(:))) / (2*numel(b1));
end