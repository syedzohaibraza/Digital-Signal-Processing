clear all;
clc;

% Specify amplitudes
a0 = -1;
a1 = 3;
a2 = 0.3;
an = 0;

% Specify frequencies [Hz]
f1 = 742.2;
f2 = 1055;

% Specify angles [rad]
phi1 = 0;
phi2 = pi/2;

% Specify sampling frequency
Fs=10000; 

% Specify number of samples
N=128;  % Number of samples (needs to be a power of 2 for FFT!!!)

% Specify Windowing (Hanning window = 1, rectangular window = 0)
hanning_window = 1;

% Specify length of FFT (NFFT-point FFT, padded with zeros if sampled
% signal has less than NFFT points and truncated if it has more
NFFT = N; % here simplified!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Simulated continuous time domain signal"
t=(0:6000)*10e-6;
x=a0+a1*cos(2*pi*f1*t+phi1)+a2*cos(2*pi*f2*t+phi2)+an*randn(size(t));

% "Simulated sampled signal"
T=1/Fs; % sampling period
n=0:N-1;
xs=a0+a1*cos(2*pi*f1*n*T+phi1)+a2*cos(2*pi*f2*n*T+phi2)+an*randn(size(n*T));

if hanning_window > 0 % windowing
  w= hann(N); % hanning window
else
  w=ones(N,1); % rectangular window
end
xsw=xs.*w';

% DFT of sampled signal
disp('Resolution frequency of DFT [Hz]')
df=Fs/NFFT
X=fftshift(fft(xsw,NFFT));  % Nfft-point FFT, padded with zeros if xsw has less than 
                            % NFFT points and truncated if it has more
Omega=(-NFFT/2:NFFT/2-1)*pi/(NFFT/2);    % DT Frequency... -pi to pi
f=(Fs/2)*(Omega/pi);   %  Convert into equivalent Hz values


% DFT of sampled signal with zero padding to reduce the grid error
NFFTZP=2^nextpow2(NFFT)
NFFTZP=16*NFFTZP
XZP=fftshift(fft(xsw,NFFTZP));
OmegaZP=(-NFFTZP/2:NFFTZP/2-1)*pi/(NFFTZP/2);    % DT Frequency... -pi to pi
fZP=(Fs/2)*(OmegaZP/pi);   %  Convert into equivalent Hz values

% Time domain and frequency domain plots
subplot(3,1,1)
plot(t,x);
grid on;
V=axis;
axis([0 N*T V(3) V(4)])
xlabel('t (sec)')
ylabel('Amplitude')

subplot(3,1,2)
stem(n,xsw,'filled')
grid on;
V=axis;
axis([0 N V(3) V(4)])
xlabel('Sample #')
ylabel('Amplitude')

subplot(3,1,3)
h=stem(f,1/sum(w)*abs(X)); % Divide by the sum of window to get the real amplitudes
%set(h,'linewidth',2)
%set(gca,'fontsize',12)
%xlabel('f  (Hz)')
%ylabel('|X_{sw}[k]|')
hold on

h=plot(fZP,1/sum(w)*abs(XZP),'r'); % Divide by the sum of window to get the real amplitudes
%set(h,'linewidth',1)
grid on;
hold off
