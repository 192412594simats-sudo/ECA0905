
clc;
clear;
close all;

%% INPUTS
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

Bc = 8;       % Coefficient word length
Bp = 16;      % Product word length
N = 200;

x = [1 zeros(1,N-1)];  % Unit impulse

%% 1. COEFFICIENT QUANTIZATION
% Signed fractional quantization
S = 2^(Bc-1);

bq = round(b*S)/S;
aq = round(a*S)/S;

% Keep denominator leading coefficient normalized
aq = aq/aq(1);

CoeffError = [b-bq, a-aq];

%% 2. IDEAL AND COEFFICIENT-QUANTIZED RESPONSE
y_ideal = filter(b,a,x);
y_coeff = filter(bq,aq,x);

e_coeff = y_ideal-y_coeff;

MaxCoeffError = max(abs(e_coeff));

if sum(e_coeff.^2) > 0
    SQNR_coeff = 10*log10( ...
        sum(y_ideal.^2)/sum(e_coeff.^2));
else
    SQNR_coeff = Inf;
end

%% 3. PRODUCT QUANTIZATION
Sp = 2^(Bp-1);

y_prod = zeros(1,N);
v = zeros(1,4);

for k = 1:N

    % Feedforward products
    p = bq*x(k);

    % Feedback products
    fb = -aq(2:5).*v;

    % Quantize products separately
    p = round(p*Sp)/Sp;
    fb = round(fb*Sp)/Sp;

    % Sum quantized products
    y_prod(k) = sum(p)+sum(fb);

    % Update filter state
    v = [y_prod(k), v(1:3)];

end

e_prod = y_coeff-y_prod;

MaxProdError = max(abs(e_prod));

if sum(e_prod.^2) > 0
    SQNR_prod = 10*log10( ...
        sum(y_coeff.^2)/sum(e_prod.^2));
else
    SQNR_prod = Inf;
end

%% 4. ZERO-INPUT LIMIT-CYCLE TEST
state = [0.01 0.01 0.01 0.01];
y_zero = zeros(1,N);

for k = 1:N

    acc = -aq(2:5)*state.';

    % Quantize output to product precision
    y_zero(k) = round(acc*Sp)/Sp;

    % Update state
    state = [y_zero(k), state(1:3)];

end

LimitAmplitude = max(abs(y_zero));
LimitSamples = find(abs(y_zero)>1e-12);

%% 5. OVERFLOW / SATURATION TEST
xo = 1.2*ones(1,N);

y_over = filter(bq,aq,xo);

% Normalized saturation to [-1,1]
y_sat = min(max(y_over,-1),1);

OverflowCount = sum(abs(y_over)>1);
MaxOverflow = max(abs(y_over));

%% 6. OUTPUT RESULTS
fprintf('\n============================================\n');
fprintf('       FIXED-POINT QUANTIZATION RESULTS\n');
fprintf('============================================\n');

fprintf('Coefficient word length = %d bits\n',Bc);
fprintf('Product word length = %d bits\n',Bp);

fprintf('\nMaximum coefficient/output error = %.6e\n', ...
    MaxCoeffError);
fprintf('Coefficient SQNR = %.2f dB\n',SQNR_coeff);

fprintf('\nMaximum product quantization error = %.6e\n', ...
    MaxProdError);
fprintf('Product SQNR = %.2f dB\n',SQNR_prod);

fprintf('\nZero-input maximum amplitude = %.6e\n', ...
    LimitAmplitude);
fprintf('Nonzero zero-input samples = %d\n', ...
    length(LimitSamples));

fprintf('\nOverflow samples = %d\n',OverflowCount);
fprintf('Maximum absolute output before saturation = %.6e\n', ...
    MaxOverflow);

%% OUTPUT PARAMETER TABLE
Test = {
    'Coefficient Quantization'
    'Product Quantization'
    'Zero-Input Test'
    'Overflow Test'
};

MaximumValue = [
    MaxCoeffError
    MaxProdError
    LimitAmplitude
    MaxOverflow
];

SQNR_dB = [
    SQNR_coeff
    SQNR_prod
    NaN
    NaN
];

T = table(Test,MaximumValue,SQNR_dB);
disp(T);

%% 7. QUANTIZED COEFFICIENTS
fprintf('\nOriginal numerator b =\n');
disp(b);

fprintf('Quantized numerator bq =\n');
disp(bq);

fprintf('Original denominator a =\n');
disp(a);

fprintf('Quantized denominator aq =\n');
disp(aq);

%% GRAPH 1: IDEAL VS COEFFICIENT-QUANTIZED RESPONSE
figure('Name','Ideal vs Quantized Response');

stem(0:N-1,y_ideal,'b','filled');
hold on;
stem(0:N-1,y_coeff,'r');

grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Ideal vs Coefficient-Quantized Response');
legend('Ideal','Quantized');

%% GRAPH 2: COEFFICIENT QUANTIZATION ERROR
figure('Name','Coefficient Quantization Noise');

plot(0:N-1,e_coeff,'LineWidth',1);

grid on;
xlabel('Sample n');
ylabel('Error');
title('Coefficient Quantization Noise');

%% GRAPH 3: PRODUCT QUANTIZATION ERROR
figure('Name','Product Quantization Noise');

stem(0:N-1,e_prod,'filled');

grid on;
xlabel('Sample n');
ylabel('Error');
title('Product Quantization Noise');

%% GRAPH 4: ZERO-INPUT LIMIT-CYCLE TEST
figure('Name','Zero-Input Limit-Cycle Test');

plot(0:N-1,y_zero,'LineWidth',1);

grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Zero-Input Limit-Cycle Test');

%% GRAPH 5: OVERFLOW / SATURATION TEST
figure('Name','Overflow and Saturation');

plot(0:N-1,y_over,'b');
hold on;
plot(0:N-1,y_sat,'r--');

grid on;
xlabel('Sample n');
ylabel('Amplitude');
title('Overflow / Saturation Test');
legend('Without Saturation','With Saturation');

%% GRAPH 6: IDEAL FILTER POLE-ZERO PLOT
figure('Name','Ideal Filter Pole-Zero Plot');

zplane(b,a);
grid on;
title('Ideal Filter Pole-Zero Plot');

%% GRAPH 7: QUANTIZED FILTER POLE-ZERO PLOT
figure('Name','Quantized Filter Pole-Zero Plot');

zplane(bq,aq);
grid on;
title('Quantized Filter Pole-Zero Plot');

%% FINAL RESULT
fprintf('\n============================================\n');
fprintf('                FINAL RESULT\n');
fprintf('============================================\n');

fprintf('1. Coefficient quantization was performed.\n');
fprintf('2. Product quantization was evaluated.\n');
fprintf('3. Zero-input behavior was tested.\n');
fprintf('4. Overflow and saturation were evaluated.\n');
fprintf('5. Pole-zero plots were generated for both filters.\n');
fprintf('============================================\n');
