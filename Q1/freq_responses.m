close all
clear


% Filter Banks

h0_6tap = [0.0352 -0.0854 -0.1350 0.4599 0.8069 0.3327];
h1_6tap = [-0.3327 0.8069 -0.4599 -0.1350 0.0854 0.0352];
g0_6tap = [0.3327    0.8069    0.4599   -0.1350   -0.0854    0.0352];
g1_6tap = [0.0352    0.0854   -0.1350   -0.4599    0.8069   -0.3327];


h0_8tap = [-0.0106    0.0329    0.0308   -0.1870   -0.0280    0.6309    0.7148    0.2304];
h1_8tap = [-0.2304    0.7148    -0.6309   -0.0280   0.1870    0.0308    -0.0329   -0.0106];
g0_8tap = [0.2304    0.7148    0.6309   -0.0280   -0.1870    0.0308    0.0329   -0.0106];
g1_8tap = [-0.0106   -0.0329    0.0308    0.1870   -0.0280   -0.6309    0.7148   -0.2304];

h0_10tap = [0.0033   -0.0126   -0.0062    0.0776   -0.0322   -0.2423    0.1384    0.7243    0.6038    0.1601];
h1_10tap = [-0.1601    0.6038    -0.7243    0.1384   0.2423   -0.0322    -0.0776   -0.0062   0.0126    0.0033];
g0_10tap = [0.1601    0.6038    0.7243    0.1384   -0.2423   -0.0322    0.0776   -0.0062   -0.0126    0.0033];
g1_10tap = [0.0033    0.0126   -0.0062   -0.0776   -0.0322    0.2423    0.1384   -0.7243    0.6038   -0.1601];

figure
hold on
grid on
x6=abs(fft(h0_6tap,1000)).^2;
y6=abs(fftshift(fft(g0_6tap,1000))).^2;

plot(x6)
plot(y6)

x8=abs(fft(h0_8tap,1000)).^2;
y8=abs(fftshift(fft(g0_8tap,1000))).^2;
plot(x8)
plot(y8)

x10=abs(fft(h0_10tap,1000)).^2;
y10=abs(fftshift(fft(g0_10tap,1000))).^2;
plot(x10)
plot(y10)

ylim([0 2])
xlim([0 500])
legend("h0 6 tap","g0 6 tap","h0 8 tap","g0 8 tap","h0 10 tap", "g0 10 tap")
title("Magnitude Square Response")
ylabel("Magnitude Squared Frequency Responses")
xlabel("Frequencies (mHz)")

figure
hold on

x6=abs(fft(h0_6tap,1000)).^2;
y6=abs(fftshift(fft(g0_6tap,1000))).^2;
plot(x6+y6);

x8=abs(fft(h0_8tap,1000)).^2;
y8=abs(fftshift(fft(g0_8tap,1000))).^2;
plot(x8+y8);

x10=abs(fft(h0_10tap,1000)).^2;
y10=abs(fftshift(fft(g0_10tap,1000))).^2;
plot(x10+y10);

legend("6 tap","8 tap","10 tap");
xlabel('Frequency');
ylabel('Magnitude Squared');
title('Perfect Reconstruction |H0(f)|^2 + |G0(f)|^2~=2');
ylim([1.95 2.05])
grid on