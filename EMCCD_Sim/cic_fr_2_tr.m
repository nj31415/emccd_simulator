%% Nicholas Jones - njones31@vt.edu
% Script to study how to translate CIC in terms of e- pixel^-1 fr^-1 to e-
% pixel^-1 tr^-1
close all;
clear;
clc;

% Study the distribution of the number of transfers in a frame per pixel
% Image size
img_w = 25; % Pixels
img_l = 25; % Pixels

num_reg = 50 + 16; % Number of register elements

num_tr = zeros(img_l, img_w);

for i = 1 : img_l
    for j = 1 : img_w
        num_tr(i, j) = i + j + num_reg;
    end
end

% Determine the mean number of transfers
mean_tr = mean(num_tr, 'all');

% Find the total number of transfers that occur in a read out
sum_tr = sum(num_tr, 'all');

% Plot a histogram of the number of transfers for all the different pixels
% in the parallel section
histogram(num_tr, 'BinWidth', 1);
xlabel('Number of Transfers');
ylabel('Count');
xline(mean_tr, 'k:', 'Label', 'Mean Transfers');