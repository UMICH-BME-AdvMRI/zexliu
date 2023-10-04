%% Adv Topics in MRI HW1 Zexuan Liu
clear; close all; clc;
%% Problem 2: Balanced and spoiled steady-state sequences
% Code borrowed and edited from Brian Hargreaves
% a) Simulate steady-state frequency response of a bSSFP sequnece

T1 = 1000; % ms
T2 = 50; % ms
TE = [2.5 5 10]; % ms
TR = [5 10 20]; % ms
flip = pi/3; % 60 degree

df = [-100:100]; % Hz

Sig = zeros(length(df),length(TE));

% TE=2.5ms, TR=5ms
for k=1:length(df)
    [Msig_2a1,Mss_2a1] = sssignal(flip,T1,T2,TE(1),TR(1),df(k));
    Sig(k,1) = Mss_2a1(1)+i*Mss_2a1(2);
end

% TE=5ms, TR=10ms
for k=1:length(df)
    [Msig_2a2,Mss_2a2] = sssignal(flip,T1,T2,TE(2),TR(2),df(k));
    Sig(k,2) = Mss_2a2(1)+i*Mss_2a2(2);
end

% TE=10ms, TR=20ms
for k=1:length(df)
    [Msig_2a3,Mss_2a3] = sssignal(flip,T1,T2,TE(3),TR(3),df(k));
    Sig(k,3) = Mss_2a3(1)+i*Mss_2a3(2);
end

figure
subplot(2,1,1);
plot(df,abs(Sig));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('TE=2.5,TR=5', 'TE=5,TR=10', 'TE=10,TR=20');

subplot(2,1,2);
plot(df,angle(Sig));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis([min(df) max(df) -pi pi]);
grid on;
legend('TE=2.5,TR=5', 'TE=5,TR=10', 'TE=10,TR=20');

%% 2-b) Modify the bSSFP sequence to generate a FLASH sequence by adding a 
% gradient spoiler along the slice selection direction

% i) Perfect spoiler: plot magnetization
T1 = 1000; % ms
T2 = 100; % ms
TE = 5; % ms
TR = 10; % ms
flip = pi/18; % 10 degree
dfreq = 0; % Assuming no off-resonance
phi = 0; % Completely dephase by end of TR
[Msig_2b1,Mss_2b1] = gssignal(flip,T1,T2,TE,TR,dfreq,phi);
sig_mag = abs(Msig_2b1);
sig_phase = angle(Msig_2b1);
% sig_mag and sig_phase is the steady state signal

% ii) Spoiler with various dephasing
phi = [2*pi 4*pi 8*pi 16*pi]; % Various dephasing cycles
sig = NaN(length(phi),2);
for k = 1:length(phi)
    [Msig_2b2,Mss_2b2] = gssignal(flip,T1,T2,TE,TR,dfreq,phi(k));
    sig(k,1) = abs(Msig_2b2);
    sig(k,2) = angle(Msig_2b2);
end
figure
subplot(2,1,1)
plot(phi/pi,sig(:,1),'o'); xlabel('Dephasing per voxel (pi)'); ylabel('Magnitude')
title('Steady state signal with spoiler gradient with various dephasing by the end of TR')
subplot(2,1,2)
plot(phi/pi,sig(:,2),'o'); xlabel('Dephasing per voxel (pi)'); ylabel('Phase')

% iii) Both gradient spoiler and RF spoiling
inc = [0:1:180]/180*pi; % phase sweep from 0-180 degree (5 degree apart)
Nex = 100; % number of excitation
Msig_2b3 = NaN(1,length(inc));
for k = 1:length(inc)
    [Msig_2b3(k),Mss_2b3] = spgrsignal(flip,T1,T2,TE,TR,dfreq,Nex,inc(k));
end
figure
plot([0:1:180],abs(Msig_2b3)); xlabel('Phase sweep (degree)'); ylabel('Mss_x_y')
% Not sure if it is correct, but there seem to be many choices to achieve
% low Mss_xy

%% 3-1&3-2) Design RF pulse and slice selective gradient to excite a slice of 5 mm
TBW = 8;
del_z = 0.005; % m
tau_rf = 0.002; % pulse duration (s)
B1_max = 25; % mT/m
slew_max = 180*1000; % T/(m*ms)
gamma = 42.58*1e6; % Hz/T
BW = TBW/(tau_rf); % bandwidth (Hz)
dT = 1/BW; % time step (s)
gradient = BW/(gamma*del_z); % gradient strength (T/m)
% build sinc pulse
xm = (9-1)/2; 
x = [-xm:xm]/xm; 
h = sinc(x*TBW/2).*(0.54+0.46*cos(pi*x)); 
h = h/sum(h);
rf = (pi/2)*h; % scaled to 90 degree
rfs = rf./(gamma*dT); % Sinc pulse
zlocs = [-2.5:0.01:2.5]/1000; % m
% df is 0 or 200 Hz
dT_blochsim = dT*1000; % ms
g = max(rfs); % mT/mm
time_ramp = gradient/slew_max; % ms

% figure plot
figure
subplot(2,1,1); plot(x,rfs*1000); title('RF waveform'); xlabel('Time (ms)'); ylabel('Amplitude (µT)')
% May still be some issue with this waveform, not exactly a sinc pulse, but
% the time step is large compared to the duration of the pulse there are
% only 9 time steps...
subplot(2,1,2); plot([0 2],gradient*ones(1,2)); title('Gradient'); xlabel('Time (ms)'); ylabel('Amplitude (mT/mm)')
% time to ramp up and down the gradient are not drawn in this plot for its
% small magnitude in time, but is calcualted above in variable time_ramp

% For 0 Hz off-resonance
df = 0; 
M_mtx_0offres = NaN(4,length(zlocs));
for k = 1:length(zlocs)
    [M_mtx_0offres(1:3,k)] = blochsim(rfs,zeros(1,length(rfs)),g*ones(1,length(rfs)),dT_blochsim,df,zlocs(k),1);
    M_mtx_0offres(4,k) = abs(M_mtx_0offres(1,k)+i*M_mtx_0offres(2,k));
end
figure
subplot(2,2,1); plot(zlocs*1000,M_mtx_0offres(1,:)); title('Mx 0Hz off res'); xlabel('Position (mm)'); ylabel('Mx')
subplot(2,2,2); plot(zlocs*1000,M_mtx_0offres(2,:)); title('My 0Hz off res'); xlabel('Position (mm)'); ylabel('My')
subplot(2,2,3); plot(zlocs*1000,M_mtx_0offres(3,:)); title('Mz 0Hz off res'); xlabel('Position (mm)'); ylabel('Mz')
subplot(2,2,4); plot(zlocs*1000,M_mtx_0offres(4,:)); title('Mxy 0Hz off res'); xlabel('Position (mm)'); ylabel('Mxy')

% For 200 Hz off-resonance
df = 200; 
M_mtx_200offres = NaN(4,length(zlocs));
for k = 1:length(zlocs)
    [M_mtx_200offres(1:3,k)] = blochsim(rfs,zeros(1,length(rfs)),g*ones(1,length(rfs)),dT_blochsim,df,zlocs(k),1);
    M_mtx_200offres(4,k) = abs(M_mtx_200offres(1,k)+i*M_mtx_200offres(2,k));
end
figure
subplot(2,2,1); plot(zlocs*1000,M_mtx_0offres(1,:)); title('Mx 200Hz off res'); xlabel('Position (mm)'); ylabel('Mx')
subplot(2,2,2); plot(zlocs*1000,M_mtx_0offres(2,:)); title('My 200Hz off res'); xlabel('Position (mm)'); ylabel('My')
subplot(2,2,3); plot(zlocs*1000,M_mtx_0offres(3,:)); title('Mz 200Hz off res'); xlabel('Position (mm)'); ylabel('Mz')
subplot(2,2,4); plot(zlocs*1000,M_mtx_0offres(4,:)); title('Mxy 200Hz off res'); xlabel('Position (mm)'); ylabel('Mxy')

%% 3-3) Different flip angle
rf_30 = (pi/6)*h; % scaled to 30 degree
rfs_30 = rf_30./(gamma*dT); % Sinc pulse
figure
subplot(3,1,1); plot(x,rfs_3_3*1000); title('RF pulse for 30 degree flip angle'); xlabel('Time (ms)'); ylabel('Amplitude (µT)')
rf_10 = (pi/18)*h; % scaled to 10 degree
rfs_10 = rf_10./(gamma*dT); % Sinc pulse
subplot(3,1,2); plot(x,rfs_10*1000); title('RF pulse for 10 degree flip angle'); xlabel('Time (ms)'); ylabel('Amplitude (µT)')
rfs_90_ft = fftshift(fft(rfs));
rfs_10_ft = fftshift(fft(rfs_10));
subplot(3,1,3); 
plot(linspace(-BW/2,BW/2,9),abs(rfs_90_ft)); title('FT of RF pulses'); xlabel('Frequency (Hz)'); ylabel('Amplitude (µT)')
hold on; plot(linspace(-BW/2,BW/2,9),abs(rfs_10_ft));
legend('90 degree flip angle','10 degree flip angle')
% Amplitude of the pulse with the smaller flip angle is smaller
% Bandwidth and function are overall the same for these flip angles

%% 3-4) Rephasing gradient
% Rephasing gradient is to correct dephasing happening in the transverse
% plane after slice select gradient.
% Still need to figure out incorporating rephasing gradient... Need to add
% relaxation to my blochsim function because previously T1 and T2 were long
% in 3-1 and 3-2...