%% Nicholas Jones - njones31@vt.edu
% Script for simulating EM-register gain in an EMCCD for various input
% photon levels. Based on Photon Counting Strategies with Low-Light-Level
% CCDs by Basden et al. 2003
% close all;
clear;
clc;

ph_in = 1 : 2;          % Photons - input photon levels
disp_mult = 2;          % Controls the x-axis length of the PDF plot
g = 1000;                % Mean gain. Input with calculation of p
r = 50;                % Number of multiplication elements
% p = 0.001148253;       % Multiplication probability ( 0.01 ~< p ~< 0.02).
                        % Input with calculation of g
p = nthroot(g, r) - 1;

rms_read_noise = 9.2;   % e- - Detector read noise at 1 MHz

% g = (1 + p)^r;          % Mean gain

x_in = ones(100000, length(ph_in)) .* ph_in;
x_out = x_in;
p_x = (1 : g * 4 * disp_mult)';
p_dist = zeros(length(p_x), min(size(x_in)));

px = @(x, n, g) x.^(n - 1) .* exp(-x ./ g) ./ (g.^n .* factorial(n - 1));

% Calculate and plot the numerical output PDF
figure();
%subplot(4, 1, 1);

for i = 1 : length(ph_in)
    parfor j = 1 : length(x_in)
        for k = 1 : r
            x_out(j, i) = x_out(j, i) + binom_rnd(x_out(j, i), p);
        end
    end

    p_dist(:, i) = px(p_x, ph_in(i), g);

%     histogram(x_out(:, i), 100, 'Normalization', 'pdf', 'FaceAlpha', ...
%         0.2, 'FaceColor', get_color(i), 'EdgeColor', 'none');
    histogram(x_out(:, i), 'BinMethod', 'integers', 'Normalization', ...
        'pdf', 'FaceAlpha', 0.2, 'FaceColor', get_color(i), ...
        'EdgeColor', 'none');
    hold on;
    plot(p_x, p_dist(:, i), '.', 'Color', get_color(i));
end

title('Output Probability Distribution for Given Electron Inputs');
xlabel('Output (e^-)');
ylabel('Probability');
xlim([0 length(p_x)]);
xline(g * ph_in, 'k--');
xline(mean(x_out, 1), 'r--');
mean_error = 100 * (mean(x_out, 1) - (g * ph_in)) ./ (g * ph_in);
leg1 = legend('1', '', '2', '', '3', '', '4', '', '', '', '', '', '', ...
    '', '', '');
title(leg1, 'Electron Input');

%% Function to generate binomial random numbers. Faster method than
% binornd from the Statistics Toolbox. Implements gpuArray for
% number of trials greater than 20000 (GPU based processing becomes
% faster at this point according to some basic comparison testing
% using timeit).
% Inputs;
% n:    Integer - number of trials
% p:    Float - probability of success
% Outputs:
% rnd:  Integer - random number from the binomial distribution
function rnd = binom_rnd(n, p)
    if n > 20000
        rnd = gather(sum(gpuArray.rand(n, 1) < p));
    else
        rnd = sum(rand(n, 1) < p);
    end
end

% % Plot model PDFs along with the RMS noise PDF
% subplot(4, 1, 2);
% plot(p_x, pdf('Normal', p_x, 0, rms_read_noise), 'Color', ...
%     get_color(i + 1));
% hold on;
% 
% for i = 1 : length(ph_in)
%     plot(p_x, p_dist(:, i), '.', 'Color', get_color(i));
% end
% 
% title('Output Probability Distribution for Given Photon Inputs');
% xlabel('Output (e^-)');
% ylabel('Probability');
% legend('Read Noise', '1', '2', '3', '4');
% xlim([0 length(p_x)]);
% ylim([0 1e-3]);
% set(gca, 'XScale', 'log');
% 
% % Plot the Poisson PDF for various light levels
% poi_x = (0 : 10)';
% poi_dist = zeros(length(poi_x), min(size(x_in)) + 1);
% 
% subplot(4, 1, 3);
% poi_dist(:, end) = poisspdf(poi_x, 0.1);
% plot(poi_x, poi_dist(:, end), '*--', 'Color', get_color(0));
% hold on;
% 
% for i = 1 : length(ph_in)
%     poi_dist(:, i) = poisspdf(poi_x, ph_in(i));
%     plot(poi_x, poi_dist(:, i), '*--', 'Color', get_color(i));
% end
% 
% title('Poisson Probability Distribution for Given Mean Photon Inputs');
% xlabel('Incident Photons (photons)');
% ylabel('Probability');
% legend('0.1', '1', '2', '3', '4');
% xlim([0 length(poi_x) - 1]);
% 
% % Simulate Poisson distributed input to multiplication register
% poi_x_in = poissrnd(0.1, 10000, 1);
% poi_x_out = poi_x_in;
% 
% subplot(4, 1, 4);
% 
% for j = 1 : length(poi_x_in)
%     for k = 1 : r
%         poi_x_out(j) = fx(poi_x_out(j), p);
%     end
% end
% 
% [~, edges] = histcounts(log10(poi_x_out));
% histogram(poi_x_out, 10.^edges, 'Normalization', 'pdf', 'FaceAlpha', ...
%     0.2, 'FaceColor', get_color(1), 'EdgeColor', 'none');
% hold on;
% plot(1 : 100, pdf('Normal', 1 : 100, 0, rms_read_noise), 'Color', ...
%      get_color(0));
% 
% title(['Output Probability Distribution for Mean Photon Input of 0.1 '...
%     'photons (Poisson Distribution)']);
% xlabel('Output (e^-)');
% ylabel('Probability');
% ylim([1e-7 inf]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% 
% % Apply thresholding and look at percent error between input and output
% % signals for photon counting
% tot_sig = sum(poi_x_in);
% read_noise_vec = rms_read_noise * randn(size(poi_x_out));
% pxo_noise = poi_x_out + read_noise_vec;
% thresh_sig = sum(pxo_noise >= (rms_read_noise * 5));
% sig_err = 100 * (thresh_sig - tot_sig) / tot_sig;

% Apply correction. thresh_sig must be in terms of flux to work (divide by
% the length of poi_x_in, which is the 'time' unit in this case)
% corr_sig = -log(1 - thresh_sig);
% sig_err_corr = 100 * (corr_sig - tot_sig) / tot_sig;

% Apply thresholding and look at percent error between input and output
% signals for multi-thresholding strategy


%% Function to get a marker color based on the supplied index.
% Inputs:
% idx   : Int, style number, will be modded to access valid index
%         in marker_color vector
% Outputs:
% m_c   : character vector, the marker color to use
function m_c = get_color(idx)
    marker_color = {'#0072BD', '#EDB120', '#77AC30', '#A2142F', ...
        '#D95319'};
    m_c = marker_color(mod(idx, length(marker_color)) + 1);
    m_c = m_c{:};
end