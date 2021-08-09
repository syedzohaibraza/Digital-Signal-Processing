%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  I. Signal Access and Exploration  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clc;
clear all;

[x,Fs]=wavread('guitar1.wav');
x=x.';   % convert into row vector

X1=fftshift(fft(x(20000+(1:16384)),65536));

f=(-32768:32767)*Fs/65536;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  II. Adding A High Frequency Interference  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=2*pi*(10000/Fs);  % convert 10 kHz into DT frequency
N=length(x);
n=0:(N-1);

x_10=x+0.2*cos(omega*n);

X_10_1=fftshift(fft(x_10(20000+(1:16384)),65536));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  III. Lowpass Filter Design  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lowpass Filter Design Specifications:
% ·	Lowpass filter 
% ·	Passband cutoff frequency = 7 kHz
% ·	Stopband cutoff frequency = 9 kHz
% ·	Sampling Frequency = 44.1 kHz  (frequencies of interest 0 to 22.05 kHz)
% ·	At least 60 dB of stopband attenuation
% ·	No more than 1 dB passband ripple

rp=1;    % specify passband ripple in dB
rs=60;    % specify stopband attenuation in dB
f_spec=[7000 9000];   % specify passband and stopband edges in Hz
AA=[1 0];    %%%  specfies that you want a lowpass filter
dev=[(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];   % compute a parameter needed by design routine
Fs=44.1e3;

[N,fo,ao,w]=firpmord(f_spec,AA,dev,Fs);     % estimates filter order and gives other design parameters

b=firpm(N,fo,ao,w);     % Computes the designed filter coefficients in vector b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  IV. Remove Interference w/ Filter  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_10_out=filter(b,1,x_10);   %%% filter the signal with the designed filter

%%%% Compare in the frequency domain
X_10_out_1=fftshift(fft(x_10_out(20000+(1:16384)),65536));

figure
subplot(3,1,1)
plot(f/1e3,20*log10(abs(X_10_1)));
title('DFT of Signal w/ Interference')
axis([-20 20 -60 90])
xlabel('f (kHz)')
ylabel('|DFT|  (dB)')
grid
subplot(3,1,2)
plot(f/1e3,20*log10(abs(X_10_out_1)));
title('DFT of Filtered Signal')
axis([-20 20 -60 90])
xlabel('f (kHz)')
ylabel('|DFT|  (dB)')
grid
subplot(3,1,3)
plot(f/1e3,20*log10(abs(X1)));
title('DFT of Original Signal')
axis([-20 20 -60 90])
xlabel('f (kHz)')
ylabel('|DFT|  (dB)')
grid

figure
t=(0:49999)*(1/Fs);
subplot(3,1,1)
plot(t,x_10(1:50000),'r')
title('Signal w/ Interference')
xlabel('t  (s)')
ylabel('signal x_10[n]')
grid
axis([0.671 0.68 -0.3 0.3])
subplot(3,1,2)
plot(t,x(1:50000),'b',t,x_10_out(1:50000),'m--')
title('Filtered Signal and Original')
xlabel('t  (s)')
ylabel('signal x_{10_out}[n] and x[n]')
grid
axis([0.671 0.68 -0.3 0.3])
legend('Original','Filtered')
subplot(3,1,3)
%%%%% Make a plot that accounts for the delay in the filtered signal
%%%%   this allows a better comparison
%%% For an odd filter order N the delay is (N-1)/2
plot(t,x(1:50000),'b',t,x_10_out(21+(1:50000)),'m--')
title('Filtered Signal (delay-corrected) and Original')
xlabel('t  (s)')
ylabel('signal x_{10_out}[n] and x[n]')
grid
axis([0.671 0.68 -0.3 0.3])
legend('Original','Filtered')

sound(x_10,Fs)
input('hit return')
sound(x_10_out,Fs)


