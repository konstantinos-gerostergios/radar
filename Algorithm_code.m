close all
clear all
clc
%% Signal creation
%parameters
T=1e-3;
fs=100e3;
DC=0.1;
Pt=1;
R0=1000;
c = 3e8;
G = 1;
L = 1;
l = 1;
w = 0;
ti = DC * T;

t = 0:1/fs:T;  % Time vector
pulse = rectpuls(t - (DC * T)/2,DC*T); %pulse transmitted
x = pulse;

stem(t, pulse);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Real part of the signal');
grid on;
ylim([-0.2 1.2]);

%% Compute DFT
K = length(pulse);  % Number of samples
k = -K/2:(K/2-1);
pulse_fft = fft(pulse,K);

figure;
stem(k/K,fftshift(abs(pulse_fft)));
xlabel('Signal frequency/Sampling frequency (f/fs)');
ylabel('DFT Amplitude');
title('Magnitude Spectrum');
grid on;
ylim([0 12]);

%% Plot the 1D data cube

% Plotting
figure;
imagesc([0 0], t, abs(pulse)'); %[0 0] in order to have 1D 
title('1D Data Cube of the Signal');
xlabel('Slow time (s)');
ylabel('Fast Time (s)');

%% Radar equation
%n = wgn(1,K,-50)

Pr = (Pt*G^2*l^2)/((64*pi^3)*R0^4*L)

pulse = zeros(size(t));
to = (2*R0)/c;    % Delay time
pulse((t >= to) & (t < to + DC*T)) = 1;
y = sqrt(Pr)*pulse; %+ n

stem(t,y);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Real part of the signal');
grid on;

figure;
imagesc([0 0], t, y');
title('1D Data Cube of the Signal');
xlabel('Slow time (s)');
ylabel('Fast Time (s)');


%% 
%% MF
% Matched Filter for PRI = 1e-3, DC = 0.1, fs = 100e3

Grx = rectpuls((DC * T)/2 - t,DC * T);
MF = conv(pulse,conj(Grx));
figure;
stem(t, MF(1:101));
title('Output of Matched Filter');
xlabel('Time (s)');
ylabel('Amplitude');

%Datacube of MF
figure;
imagesc([0 0], t, MF(1:101)');
title('1D Data Cube of the Signal');
xlabel('Slow time (s)');
ylabel('Fast Time (s)');
colorbar;

%% RCS

y = rcs_model(x, R0, Pt,G, l ,L, c, t, ti, w);

% Plot the results
figure
subplot(2,1,1);
stem(t, x);
title('Transmit Pulse');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
stem(t, y);
title('Received Signal with RCS Model');
xlabel('Time (s)');
ylabel('Amplitude');


%% 

% Calculate range resolution (∆R)
delta_R = (c * DC * T) / 2

% Calculate bin resolution
bin_resolution = 1/ fs

% Calculate the number of range bins (L)
L = T * fs

%% Algorithm to find range bin
figure;
imagesc([0 0], t, y');
xlabel('Time (s)');
ylabel('Range Bin');
title('Data Cube');
colorbar;

t0  = (2 * R0) / c;
ti = t0*fs;
upper_lim = ceil(t0 * fs);
lower_lim = upper_lim - 1;
% Display the results
disp(['Possible range of received signal: ' num2str(lower_lim) '/fs to ' num2str(upper_lim) '/fs']);
disp(['The target is in the '  num2str(upper_lim) ' range bin'])


%% range bin for R0 = 15km
new_R0 = 15000;

new_y = rcs_model(x, new_R0, Pt,G, l ,L, c, t, ti, w);

figure;
stem(t, new_y);
title('Received Signal with RCS Model');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
imagesc([0 0], t, new_y');
xlabel('Time (s)');
ylabel('Range Bin');
title('Data Cube');
colorbar;

new_t0  = (2 * new_R0) / c;
new_ti = new_t0*fs;
upper_lim = ceil(new_t0 * fs);
lower_lim = upper_lim - 1;
% Display the results
disp(['Possible range of received signal: ' num2str(lower_lim) '/fs to ' num2str(upper_lim) '/fs']);
disp(['The target is in the '  num2str(upper_lim) ' range bin'])

%% ΔR/2
% ΔR = (c * τ)/2,  ΔR' = (c * τ)/4
% We have to divide ti by 2 so DC is being divided by 2
ti = (DC * T )/2;

pulse = rectpuls(t - ti/2,(DC*T)/2);  % Generate pulse based on condition
x = pulse;

y = rcs_model(x, R0, Pt,G, l ,L, c, t, ti, w);

% Plot the results
figure
subplot(2,1,1);
stem(t, x);
title('Transmit Pulse');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
stem(t, y);
title('Received Signal with RCS Model');
xlabel('Time (s)');
ylabel('Amplitude');

%% 
%% Doppler and velocity resolution
%%parameters
N = T * fs;
M = 20;
fc = 1e9;

dfd = 1 / (M * T)
dfv = (c * dfd) / (2 * fc)

%% max doppler, range and velocity
PRF = 1/T;


% I will have the maximum FD when FD is equal to PRF/2
% From PRF/2 and beyond i will have aliasing

max_fd = PRF / 2 
max_fv = (c * max_fd) / (2 * fc)
maxR = ((T- (DC * T)) * c) / 2

%% Ερώτημα γ DFT

%y2D = zeros(M, fs*T);
t0  = (2 * R0) / c;
ntarget = ceil(t0 * fs)
ti  = DC * T;
y2D = zeros(M, fs*T);

for m = 1:M
    for n = 1:(fs*T)
        y2D(m,n) =  sqrt(Pr)*rectpuls(((n-1)/fs)-(2*R0/c) -(ti/2),ti)*exp((-1i*4*pi*fc*R0*fc)/c)*exp(1i*2*pi*v*fc*(((n-1)/fs)+(m-1)*T)/c); %+ wgn(1,1,10*log10(0.0000000000000000000000045)) ;
    end
end

figure;
imagesc(abs(y2D)'); 
title('2D Datacube');
xlabel('slow time (s)');
ylabel('fast time (s)');

figure;
y2D_2 = fft(y2D(:,ntarget+1),30);
stem(fftshift(abs(y2D_2)));
title('Magnitude Spectrum');
xlabel('Signal Frequency/Sampling Frequency (f/fs)');
ylabel('DFT Amplitude');


%% Speed and doppler bin
fd = (2 * v * fc) / c;
doppler_bin = round(fd/dfd)

%% Range-Doppler Profile 

real_y2D = zeros(M,fs*T);
for n = 1:(fs*T)
    real_y2D(:,n) = fft(y2D(:,n));
end

figure;
imagesc(fftshift(abs(real_y2D))');
title('Range-Doppler Profile')
xlabel('velocity (m/s)')
ylabel('Range (m)')


%% SNR for 3D 
%% 
L = 10;
d = 0.25;
fi = pi/3; % A random angle
t0  = (2 * R0) / c;
ntarget = ceil(t0 * fs);
ti  = DC * T;
y3D = zeros(fs*T, M ,L);
SNR_dB = 15; % Set desired SNR in dB
SNR = 10.^(SNR_dB / 10); % Convert SNR from dB to linear scale
N = fs*T;

for l = 1:L
    for m = 1:M
        for n = 1:(fs*T)
            noise = sqrt(Pr/SNR)*randn(1,1); %
            y3D(n,m,l) = noise + sqrt(Pr)*rectpuls(((n-1)/fs)-(2*R0/c) -(ti/2),ti)*exp((-1i*4*pi*fc*R0*fc)/c)*exp(1i*2*pi*v*fc*(((n-1)/fs)+(m-1)*T)/c)*exp((-1i*2*pi*d*(l - 1)*cos(fi)*fc)/c);
        end
    end
end

[N, M, L] = size(y3D);
[m_slice, n_slice, l_slice] = meshgrid(1:M, 1:N, 1:L);

disp(size(n_slice));
disp(size(m_slice));
disp(size(l_slice));

% Plot using slice
figure;
slice(m_slice, n_slice, l_slice, abs(y3D), 1:M, 1:N, 1:L);
xlabel('m');  
ylabel('n');
zlabel('l');
title('3D imagesc Plot');
colorbar;

%% Azimuth Cut
%
 
theta = 0; %change theta to change the angle
for angle = 1:181
    angle_rad = (angle - 1 - theta) * pi / 180; % Convert degrees to radians
    f1_dot = cos(angle_rad) * fc / c;
    azimuth(angle) = sin(pi * L * d * f1_dot) / sin(pi * d * f1_dot);
end

% Plot
azimuth_angle = mag2db(abs(azimuth)/max(abs(azimuth)));
plot(azimuth_angle)
xlabel('Azimuth Angle (degrees)')
ylabel('Normalized Power (db)')
xlim([0 180]);
title('Azimuth Cut')


function y = rcs_model(x, R0, Pt, G, l, L, c, t, ti, w)
    Pr = (Pt*G^2*l^2)/((64*pi^3)*R0^4*L)

    % Calculate delay based on the distance
    delay = (2 * R0) / c;

    num_samples_to_shift = ceil(delay / (t(2) - t(1))) ;

    % Apply the delay by shifting the pulse vector
    y = circshift(x, num_samples_to_shift);

    % Apply RCS model
    y = sqrt(Pr) * y;
    % Add noise
    y = y + w;
end






