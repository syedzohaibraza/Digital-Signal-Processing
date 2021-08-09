function DFT_Relations(Fs,N,N_zp)
%%% A routine showing CTFT, DTFT_inf, DTFT_N, and DFT Relations
%
% Inputs:  Fs = sampling rate in Hz
%          N = number of collected samples
%          N_zp = DFT size after zero-padding
%  
%  Output:  Plots of CTFT, DTFT_inf, DTFT_N, and DFT



%%% Set parameters if user doesn't specify them
if nargin==0
    Fs=30e3;
    N=8;
    N_zp=8;
    %%% these values show significant amounts of all three error types:
            %  Aliasing Error (see 2nd plot from top)... Increase Fs to reduce
            %  Leakage Error (see 3rd plot from top)... Increase N to reduce
            %  Grid Error:  (see bottom plot)... Increase N_zp (and/or N) to reduce 
end

%%%% Compute the Theoretical CTFT
b=2*pi*1000;
f=-200000:100:200000;
CTFT=1./(j*2*pi*f + b);   %%% from table for x(t) = exp(-bt)u(t)
subplot(4,1,1)
h=plot(f/1e3,abs(CTFT),'r--');
set(h,'linewidth',2);
set(gca,'fontsize',12);
xlabel('f   (kHz)')
ylabel('|CTFT(f)|')
axis_vals=axis;
axis_vals(1:2)=[-3 3]*(Fs/2)/1000;
axis(axis_vals);

%%%%  Compute the Theoretical DTFT_inf
T=1/Fs;   
a=exp(-b*T);   %%% computes exponential decay rate of sampled signal
omega=-3*pi:0.01:3*pi;
DTFT_inf=1./(1 - a*exp(-j*omega));   %%% from table for x[n] = a^n u[n]
subplot(4,1,2)
h=plot(omega/pi,abs(T*DTFT_inf));
set(h,'linewidth',2);
set(gca,'fontsize',12);
xlabel('\Omega/\pi  (Normalized rad/sample)')
ylabel('|T*DTFT_{inf}(\Omega)|')
hold on
h=plot(f/(Fs/2),abs(CTFT),'r--');
set(h,'linewidth',2);
axis_x([-3 3])
hold off

%%%%  Compute the Theoretical DTFT_N
DTFT_N=(1-(a*exp(-j*omega)).^N)./(1 - a*exp(-j*omega));    
subplot(4,1,3)
h=plot(omega/pi,abs(T*DTFT_inf));
set(h,'linewidth',2);
set(gca,'fontsize',12);
hold on
h=plot(omega/pi,abs(T*DTFT_N),'m');
set(h,'linewidth',2);
hold off
xlabel('\Omega/\pi  (Normalized rad/sample)')
ylabel('|T*DTFT_{N}(\Omega)|')


%%%% Now do the DFT processing that could be done in practice on the N samples
n=0:(N-1);
nT=n*T;
x_N=exp(-b*nT);
%% Next command: 
   %% fft computes fft with zero-padding to N_zp
   %% fftshift moves the points for [pi,2*pi) down to [-pi, 0) so we cover [-pi,pi)
DFT_N=fftshift(fft(x_N,N_zp));   
%% now compute the frequencies where DFT points are
omega_k= ((-N_zp/2):((N_zp/2)-1))*pi/(N_zp/2);  
subplot(4,1,4)
h=plot(omega/pi,abs(T*DTFT_N),'m');
set(h,'linewidth',2);
set(gca,'fontsize',12);
hold on
h=plot(omega_k/pi,abs(T*DFT_N),'bo');
set(h,'linewidth',2);
hold off
xlabel('\Omega/\pi  (Normalized rad/sample)')
ylabel('|T*DFT_{N}(\Omega)|')
