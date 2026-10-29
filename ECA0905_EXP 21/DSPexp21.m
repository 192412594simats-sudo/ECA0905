clc;
clear;
close all;

%% LMS AND NLMS: ECHO CANCELLATION AND EQUALIZATION

rng(1);

N = 5000;                  % Number of samples
filterLength = 32;         % Adaptive filter length
muLMS = 0.01;              % LMS step size
muNLMS = 0.5;              % NLMS step size
epsilon = 1e-6;             % Prevent division by zero
noisePower = 1e-4;

%% 1. ADAPTIVE ECHO CANCELLATION
% Far-end signal and unknown echo path
farEnd = randn(N, 1);
echoPath = [0.8; 0.5; -0.3; 0.2; zeros(filterLength-4, 1)];
echo = filter(echoPath, 1, farEnd);
nearEnd = 0.1 * randn(N, 1);
microphone = echo + nearEnd;

% LMS and NLMS adaptive filters
wLMS = zeros(filterLength, 1);
wNLMS = zeros(filterLength, 1);
outLMS = zeros(N, 1);
outNLMS = zeros(N, 1);
errLMS = zeros(N, 1);
errNLMS = zeros(N, 1);

for n = filterLength:N
    u = farEnd(n:-1:n-filterLength+1);

    % LMS
    outLMS(n) = wLMS' * u;
    errLMS(n) = microphone(n) - outLMS(n);
    wLMS = wLMS + muLMS * u * errLMS(n);

    % NLMS
    outNLMS(n) = wNLMS' * u;
    errNLMS(n) = microphone(n) - outNLMS(n);
    wNLMS = wNLMS + ...
        (muNLMS / (epsilon + u' * u)) * u * errNLMS(n);
end

%% 2. EQUALIZATION
% Transmit signal through an unknown channel
tx = randn(N, 1);
channel = [1; 0.5; -0.25; zeros(filterLength-3, 1)];
rx = filter(channel, 1, tx) + sqrt(noisePower) * randn(N, 1);

% Desired equalizer output (delay allows causal equalization)
delay = 3;
desired = [zeros(delay, 1); tx(1:N-delay)];

wEqLMS = zeros(filterLength, 1);
wEqNLMS = zeros(filterLength, 1);
outEqLMS = zeros(N, 1);
outEqNLMS = zeros(N, 1);
errEqLMS = zeros(N, 1);
errEqNLMS = zeros(N, 1);

for n = filterLength:N
    u = rx(n:-1:n-filterLength+1);

    % LMS equalizer
    outEqLMS(n) = wEqLMS' * u;
    errEqLMS(n) = desired(n) - outEqLMS(n);
    wEqLMS = wEqLMS + muLMS * u * errEqLMS(n);

    % NLMS equalizer
    outEqNLMS(n) = wEqNLMS' * u;
    errEqNLMS(n) = desired(n) - outEqNLMS(n);
    wEqNLMS = wEqNLMS + ...
        (muNLMS / (epsilon + u' * u)) * u * errEqNLMS(n);
end

%% 3. PERFORMANCE RESULTS
smoothLength = 100;
mseEchoLMS = movmean(errLMS.^2, smoothLength);
mseEchoNLMS = movmean(errNLMS.^2, smoothLength);
mseEqLMS = movmean(errEqLMS.^2, smoothLength);
mseEqNLMS = movmean(errEqNLMS.^2, smoothLength);

evalIdx = floor(N/2):N;

fprintf('\n===== LMS AND NLMS PERFORMANCE =====\n');
fprintf('Echo cancellation final MSE:\n');
fprintf('  LMS  : %.6g\n', mean(errLMS(evalIdx).^2));
fprintf('  NLMS : %.6g\n', mean(errNLMS(evalIdx).^2));
fprintf('Equalization final MSE:\n');
fprintf('  LMS  : %.6g\n', mean(errEqLMS(evalIdx).^2));
fprintf('  NLMS : %.6g\n', mean(errEqNLMS(evalIdx).^2));

%% 4. PLOT CONVERGENCE
figure('Name', 'Adaptive Algorithm Performance');

subplot(2,1,1);
plot(10*log10(mseEchoLMS + eps), 'b');
hold on;
plot(10*log10(mseEchoNLMS + eps), 'r');
grid on;
title('Echo Cancellation: Learning Curves');
xlabel('Sample');
ylabel('Smoothed Error Power (dB)');
legend('LMS', 'NLMS');

subplot(2,1,2);
plot(10*log10(mseEqLMS + eps), 'b');
hold on;
plot(10*log10(mseEqNLMS + eps), 'r');
grid on;
title('Equalization: Learning Curves');
xlabel('Sample');
ylabel('Smoothed Error Power (dB)');
legend('LMS', 'NLMS');

%% 5. PLOT ECHO CANCELLATION OUTPUT
figure('Name', 'Echo Cancellation Results');

plot(1:500, microphone(1:500), 'k');
hold on;
plot(1:500, errLMS(1:500), 'b');
plot(1:500, errNLMS(1:500), 'r');
grid on;
title('Microphone Signal and Residual Echo');
xlabel('Sample');
ylabel('Amplitude');
legend('Microphone Signal', 'LMS Residual', 'NLMS Residual');