
clc;
clear;
close all;

%% FIR FILTER PARAMETERS

b = [0.0625 0.2500 0.3750 0.2500 0.0625];
Fs = 1000;
N = 128;

x = [1 zeros(1,N-1)];
n = 0:N-1;

%% 1. DIRECT REALIZATION

yd = filter(b,1,x);

%% 2. CASCADE REALIZATION

z = roots(b);

yc = x;

for k = 1:length(z)
    yc = filter([1 -z(k)],1,yc);
end

% Apply overall gain
yc = b(1)*yc;

% Remove negligible imaginary components
yc = real(yc);

%% 3. SYMMETRY AND LINEAR PHASE

symErr = max(abs(b-fliplr(b)));
D = (length(b)-1)/2;

if symErr < 1e-12
    symmetry = 'Symmetric';
else
    symmetry = 'Not Symmetric';
end

%% 4. FREQUENCY, PHASE AND GROUP DELAY

[H,w] = freqz(b,1,1024);
phase = unwrap(angle(H));
gd = grpdelay(b,1,1024);

%% 5. NUMERICAL COMPARISON

err = max(abs(yd-yc));

peak = max(abs(yd));

if peak > 0
    dev = (err/peak)*100;
else
    dev = 0;
end

%% DISPLAY RESULTS

fprintf('\n============================================\n');
fprintf('       FIR FILTER REALIZATION RESULTS\n');
fprintf('============================================\n');

fprintf('\nFIR Coefficients:\n');
disp(b);

fprintf('Sampling Frequency = %d Hz\n',Fs);
fprintf('Number of Samples = %d\n',N);

fprintf('Coefficient Symmetry = %s\n',symmetry);
fprintf('Symmetry Error = %.3e\n',symErr);
fprintf('Group Delay = %.2f samples\n',D);

fprintf('\nMaximum Error = %.3e\n',err);
fprintf('Percentage Deviation = %.3e %%\n',dev);

%% FIRST 10 IMPULSE RESPONSE VALUES

fprintf('\n============================================\n');
fprintf('     FIRST 10 IMPULSE RESPONSE VALUES\n');
fprintf('============================================\n');

fprintf('%5s %15s %15s\n','n','Direct','Cascade');

for k = 1:10
    fprintf('%5d %15.8f %15.8f\n', ...
        n(k),yd(k),yc(k));
end

%% CASCADE REALIZATION ZEROS

fprintf('\n============================================\n');
fprintf('       CASCADE REALIZATION ZEROS\n');
fprintf('============================================\n');

disp(z);

%% OUTPUT PARAMETER TABLE

Parameter = {
    'Direct realization output'
    'Cascade realization output'
    'Impulse response'
    'Phase response'
    'Comparison of realizations'
};

Observation = {
    'Direct FIR filtering'
    'Cascade of first-order sections'
    'Equal to FIR coefficient sequence'
    'Linear phase'
    'Equivalent transfer-function behavior'
};

Value = {
    sprintf('[%.4f %.4f %.4f %.4f %.4f]',yd(1:5))
    sprintf('[%.4f %.4f %.4f %.4f %.4f]',yc(1:5))
    sprintf('[%.4f %.4f %.4f %.4f %.4f]',yd(1:5))
    sprintf('Group delay = %.2f samples',D)
    sprintf('Max error = %.3e; Deviation = %.3e %%',err,dev)
};

Unit_Note = {
    'Amplitude'
    'Amplitude'
    'h[n]'
    'Samples'
    'Numerical precision'
};

Output_Table = table(Parameter,Observation,Value,Unit_Note);

fprintf('\n============================================\n');
fprintf('          OUTPUT PARAMETER TABLE\n');
fprintf('============================================\n');

disp(Output_Table);

%% GRAPH 1: DIRECT VS CASCADE IMPULSE RESPONSE

figure('Name','Direct vs Cascade Impulse Response');

stem(n,yd,'b','filled');
hold on;
stem(n,yc,'r--');

grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Direct and Cascade FIR Impulse Responses');

legend('Direct Realization','Cascade Realization');

%% GRAPH 2: MAGNITUDE RESPONSE

figure('Name','FIR Magnitude Response');

plot(w/pi,abs(H),'LineWidth',1.5);

grid on;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('|H(e^{j\omega})|');
title('FIR Magnitude Response');

%% GRAPH 3: PHASE RESPONSE

figure('Name','FIR Phase Response');

plot(w/pi,phase,'LineWidth',1.5);

grid on;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (radians)');
title('Unwrapped Phase Response');

%% GRAPH 4: GROUP DELAY

figure('Name','FIR Group Delay');

plot(w/pi,gd,'LineWidth',1.5);

grid on;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Group Delay (samples)');
title('Group Delay');

%% GRAPH 5: NUMERICAL ERROR

figure('Name','Direct vs Cascade Numerical Error');

plot(n,abs(yd-yc),'LineWidth',1.5);

grid on;
xlabel('Sample n');
ylabel('Absolute Error');
title('Direct vs Cascade Numerical Error');

%% FINAL RESULT

fprintf('\n============================================\n');
fprintf('                FINAL RESULT\n');
fprintf('============================================\n');

fprintf('1. FIR filter realized using Direct structure.\n');
fprintf('2. FIR filter realized using Cascade structure.\n');
fprintf('3. Coefficients are %s.\n',symmetry);
fprintf('4. Linear-phase group delay = %.2f samples.\n',D);
fprintf('5. Maximum realization error = %.3e.\n',err);
fprintf('6. Percentage deviation = %.3e %%.\n',dev);
fprintf('7. Direct and Cascade responses are equivalent.\n');
fprintf('============================================\n');
