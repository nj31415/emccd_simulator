%% Nicholas Jones - njones31@vt.edu
% Script for simulating EM-register gain in an EMCCD for various input
% photon levels. Based on Photon Counting Strategies with Low-Light-Level
% CCDs by Basden et al. 2003
close all;
clear;
clc;

ph_in = 1;              % Photons - input photon levels
disp_mult = 2;          % Controls the x-axis length of the PDF plot
g = 1000;               % Mean gain. Input with calculation of p
r = [604 250 100 50];   % Number of multiplication elements
p = nthroot(g, r) - 1;

x_out = ones(1000000, length(r));
p_x = (1 : g * 4 * disp_mult)';
p_dist = zeros(length(p_x), min(size(x_out)));

px = @(x, n, g) x.^(n - 1) .* exp(-x ./ g) ./ (g.^n .* factorial(n - 1));

% Calculate and plot the numerical output PDF
figure();
for i = 1 : length(r)
    r_i = r(i);
    p_i = p(i);
    parfor j = 1 : length(x_out)
        for k = 1 : r_i
            x_out(j, i) = x_out(j, i) + binom_rnd(x_out(j, i), p_i);
        end
    end

    histogram(x_out(:, i), 'BinMethod', 'integers', 'Normalization', ...
        'pdf', 'FaceAlpha', 0.4, 'FaceColor', get_color(i), ...
        'EdgeColor', 'none');
    hold on;
end

p_dist(:, i) = px(p_x, ph_in, g);
plot(p_x, p_dist(:, i), 'k.');

title('Output Probability Distribution for Given Electron Inputs');
xlabel('Output (e^-)');
ylabel('Probability');
xlim([0 length(p_x)]);
xline(g * ph_in, 'k--');
xline(mean(x_out, 1), 'r--');
mean_error = 100 * (mean(x_out, 1) - (g * ph_in)) ./ (g * ph_in);
leg1 = legend(num2str(r(1)), num2str(r(2)), num2str(r(3)), num2str(r(4)));
title(leg1, 'Num Mult');

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