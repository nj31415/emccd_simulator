%% Nicholas Jones - njones31@vt.edu
% Script for generating explanatory image of the clocking process
close all;
clear;
clc;

n_img = 10;

x_grid_px = [1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13];
y_grid_integrate = [1 1 0.1 0.1 1 1 0.1 0.1];
y_grid_combine = [1 1 0.1 0.1 -0.1 -0.1 -1 -1 ];
y_grid_transition = [0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5 -0.5];
y_grid_transfer = [0 0 -1 -1 1 1 0.1 0.1];

%% Integrate
fig = figure(1);
hold on;
rectangle('Position', [2 0.1 1 0.4], 'FaceColor', 'r');
rectangle('Position', [4 0.1 1 0.15], 'FaceColor', 'r');
rectangle('Position', [6 0.1 1 0.05], 'FaceColor', 'r');
rectangle('Position', [8 0.1 1 0.3], 'FaceColor', 'r');
rectangle('Position', [10 0.1 1 0], 'FaceColor', 'r');
rectangle('Position', [12 0.1 1 0.6], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_integrate y_grid_integrate y_grid_integrate], ...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Integrate');
drawnow;
frame = getframe(fig);
im{1} = frame2im(frame);

%% Combine
fig = figure(2);
hold on;
rectangle('Position', [4 -1 1 0.4 + 0.15], 'FaceColor', 'r');
rectangle('Position', [8 -1 1 0.05 + 0.3], 'FaceColor', 'r');
rectangle('Position', [12 -1 1 0.6], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_combine y_grid_combine y_grid_combine], ...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Combine');
drawnow;
frame = getframe(fig);
im{2} = frame2im(frame);

%% Transition
fig = figure(3);
hold on;
rectangle('Position', [4 -0.5 1 0.4 + 0.15], 'FaceColor', 'r');
rectangle('Position', [8 -0.5 1 0.05 + 0.3], 'FaceColor', 'r');
rectangle('Position', [12 -0.5 1 0.6], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_transition y_grid_transition y_grid_transition],...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Transition');
drawnow;
frame = getframe(fig);
im{3} = frame2im(frame);

%% Transfer
fig = figure(4);
hold on;
rectangle('Position', [6 -1 1 0.4 + 0.15], 'FaceColor', 'r');
rectangle('Position', [10 -1 1 0.05 + 0.3], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_transfer y_grid_transfer y_grid_transfer], ...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Transfer');
drawnow;
frame = getframe(fig);
im{4} = frame2im(frame);

%% Transition
fig = figure(5);
hold on;
rectangle('Position', [6 -0.5 1 0.4 + 0.15], 'FaceColor', 'r');
rectangle('Position', [10 -0.5 1 0.05 + 0.3], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_transition y_grid_transition y_grid_transition],...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Transition');
drawnow;
frame = getframe(fig);
im{5} = frame2im(frame);

%% Combine
fig = figure(6);
hold on;
rectangle('Position', [8 -1 1 0.4 + 0.15], 'FaceColor', 'r');
rectangle('Position', [12 -1 1 0.05 + 0.3], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_combine y_grid_combine y_grid_combine], ...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Combine');
drawnow;
frame = getframe(fig);
im{6} = frame2im(frame);


%% Transition
fig = figure(7);
hold on;
rectangle('Position', [8 -0.5 1 0.4 + 0.15], 'FaceColor', 'r');
rectangle('Position', [12 -0.5 1 0.05 + 0.3], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_transition y_grid_transition y_grid_transition],...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Transition');
drawnow;
frame = getframe(fig);
im{7} = frame2im(frame);

%% Transfer
fig = figure(8);
hold on;
rectangle('Position', [10 -1 1 0.4 + 0.15], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_transfer y_grid_transfer y_grid_transfer],...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Transfer');
drawnow;
frame = getframe(fig);
im{8} = frame2im(frame);

%% Transition
fig = figure(9);
hold on;
rectangle('Position', [10 -0.5 1 0.4 + 0.15], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_transition y_grid_transition y_grid_transition],...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Transition');
drawnow;
frame = getframe(fig);
im{9} = frame2im(frame);

%% Combine
fig = figure(10);
hold on;
rectangle('Position', [12 -1 1 0.4 + 0.15], 'FaceColor', 'r');
plot(x_grid_px, [y_grid_combine y_grid_combine y_grid_combine], ...
    'LineWidth', 3);
ylim([-1.5 1.5]);
xlim([1 13]);
xline(1.25, 'k:', 'Label', 'Pixel 1');
xline(5.25, 'k:', 'Label', 'Pixel 2');
xline(9.25, 'k:', 'Label', 'Pixel 3');
yline(0, 'b:');
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
title('Combine');
drawnow;
frame = getframe(fig);
im{10} = frame2im(frame);

%% Write the output GIF
filename = 'clock_out.gif';

for idx = 1 : n_img
    [A, map] = rgb2ind(im{idx}, 256);
    if idx == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 2);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', ...
            'DelayTime', 1.5);
    end
end