clc; clear; close all;

%% IIR Filter Coefficients
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

N = 50;
x = [1 zeros(1,N-1)]; % Unit impulse
n = 0:N-1;

%% 1. Direct Form-I and Direct Form-II
y1 = filter(b,a,x); % Direct Form-I
y2 = filter(b,a,x); % Direct Form-II

%% 2. Cascade Realization
[sos,g] = tf2sos(b,a);
yc = g*x;

for k = 1:size(sos,1)
    yc = filter(sos(k,1:3),sos(k,4:6),yc);
end

%% 3. Parallel Realization
[r,p,k] = residuez(b,a);
yp = filter(k,1,x);

for i = 1:length(r)
    yp = yp + filter(r(i),[1 -p(i)],x);
end

%% 4. Display Coefficients
disp('Numerator coefficients b ='); disp(b)
disp('Denominator coefficients a ='); disp(a)

disp('Second Order Sections (Cascade) ='); disp(sos)
disp('Cascade Gain ='); disp(g)

disp('Residues (Parallel) ='); disp(r)
disp('Poles (Parallel) ='); disp(p)
disp('Direct term ='); disp(k)

%% 5. Error and Percentage Deviation
E = [max(abs(y1-y2)); ...
    max(abs(y1-yc)); ...
    max(abs(y1-yp))];

D = E/max(abs(y1))*100;

T = table(["DF-I";"DF-II";"Cascade";"Parallel"], ...
    [0;E], ...
    [0;D], ...
    'VariableNames', ...
    {'Structure','Maximum_Error','Deviation_percent'});

disp('Comparison Table:')
disp(T)

%% 6. Impulse Response Comparison
figure(1);
stem(n,y1,'filled'); hold on;
stem(n,y2,'--');
stem(n,yc,':','LineWidth',1.5);
stem(n,yp,'-.');
grid on;
xlabel('n');
ylabel('Amplitude');
title('Impulse Response of IIR Realizations');
legend('Direct Form-I','Direct Form-II', ...
    'Cascade','Parallel');

%% 7. Error Plot
figure(2);
plot(n,abs(y1-y2),'--','LineWidth',1.2); hold on;
plot(n,abs(y1-yc),':','LineWidth',1.5);
plot(n,abs(y1-yp),'-.','LineWidth',1.2);
grid on;
xlabel('n');
ylabel('Absolute Error');
title('Numerical Error');
legend('DF-I vs DF-II','DF-I vs Cascade', ...
    'DF-I vs Parallel');

%% 8. Frequency Response
figure(3);
freqz(b,a,1024);
title('Frequency Response of IIR Filter');

%% 9. Pole-Zero Plot
figure(4);
zplane(b,a);
title('Pole-Zero Plot');

%% 10. Additional Impulse Response Graph
figure(5);
stem(n,y1,'filled');
grid on;
xlabel('n');
ylabel('Amplitude');
title('Impulse Response of IIR Filter');