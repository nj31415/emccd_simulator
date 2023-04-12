%% Nicholas Jones - njones31@vt.edu
% Script for investigating the impact of vertical and horizontal clocks on
% read-out time for the detector. Based on parameters for a CCD201-20.
% Analysis assumes that the active imaging area is located at the region
% closest to the EM read out amplifier
clear;
clc;
close all;

% Detector parameters
det_len = 2048;     % Pix - detector length
det_wid = 1024;     % Pix - detector width
N = 0;              % Pix - length of the EM register overscan, corner, and
                    % gain elements

img_len = 100;      % Pix - length of the active imaging region
img_wid = 100;      % Pix - width of the active imaging region

vert_freq = 0.9e6;  % Hz - Vertical transfer frequency
horz_freq = 1e6;    % Hz - Horizontal transfer frequency

%% Study impact of horizontal frequency on detector read out
% time
horz_freq_1 = linspace(1e6, 1e7, 100);
time_1 = det_read(N, img_len, img_wid, vert_freq, horz_freq_1);

figure();
plot(horz_freq_1 * (1e-6), time_1 * (1e3), 'k*');
yline((img_len / vert_freq) * (1e3), 'k:', 'Label', ...
    'Vertical Transfer Time');
xlabel('Horizontal Frequency (MHz)');
ylabel('Read-Out Time (ms)');
ylim([0 inf]);
title('Horizontal Frequency vs Read-Out Time');

figure();
plot(horz_freq_1 * (1e-6), 100 * (img_len ./ vert_freq) ./ time_1, 'k*');
xlabel('Horizontal Frequency (MHz)');
ylabel('Percentage of Read-Out Time (%)');
ylim([0 inf]);
title(['Horizontal Frequency vs Percentage of Read-Out Performing ' ...
    'Vertical Transfer']);

%% Function to calculate the time read out the detector
function time = det_read(N, img_len, img_wid, vert_freq, horz_freq)
time = (img_len ./ vert_freq) + ((img_len .* img_wid + N) ./ horz_freq);
end