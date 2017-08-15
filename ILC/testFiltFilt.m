%% Testing filtfilt from filtering toolbox

clc; clear; close all; rng(1);
N = 50;
t = 0.01 * linspace(1,N,N);
signal = sin(2*pi*5*t) + 0.2*cos(2*pi*6*t);
noise = 0.1 * randn(1,N);

x = signal + noise;
w = 6; % cutoff freq
w_perc = w/floor(N/2);

% filter, reverse data, filter again, and reverse data again
y1 = filtfilt2(x,w_perc);

if exist('filtfilt')
    [B,A] = butter(2,w_perc);
    y2 = filtfilt(B,A,x);
    plot(t,x,t,y1,t,y2);
    legend('signal','butterworth twice','filtfilt');
else
    plot(t,x,'-',t,y1,'--');
    legend('signal','butterworth twice');
end

rms_noise = norm(noise,2)
rms_filt = norm(y1 - signal,2)