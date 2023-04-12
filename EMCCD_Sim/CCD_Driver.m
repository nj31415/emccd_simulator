%% Nicholas Jones - njones31@vt.edu
% Driver script for the detector.
close all;
clear;
clc;

% Utilize an example iQE curve
wl_ref_qe = 1 : 0.01 : 1100;
iqe_ref = sin(2 * pi * wl_ref_qe / 1500 - pi / 2) + ...
    0.5 * sin(4 * pi * wl_ref_qe / 1000) + ...
    0.01 * sin(20 * pi * wl_ref_qe / 500).^2;
% iqe_ref = ones(size(wl_ref_qe));
iqe_ref = abs(iqe_ref) / max(iqe_ref);

% Choose between EM (1) and standard read out (0)
read_opt = 0;

figure();
plot(wl_ref_qe, iqe_ref);
xlim([wl_ref_qe(1) wl_ref_qe(end)]);
ylim([0 1]);
xlabel('Wavelegnth, nm');
ylabel('iQE');

% Instantiate the detector
par_len = 12;
par_wid = 13;
num_std_os = 2;
num_ad_elem = 2;
num_mult = 10;
em_read_noise = 30;
em_read_bias = 150;
std_read_noise = 3;
std_read_bias = 15;

test_det = Detector(par_wid, par_len, num_std_os, num_ad_elem, ...
    num_mult, em_read_noise, em_read_bias, std_read_noise, std_read_bias);

% Define pixel properties
fwc_par = 80000;
fwc_ser = 360000;
cte = 0.999999;
dcr = 0.0015;
cicr_par = 1e-5;
cicr_srl = 1e-4;

% Calculate multiplication probability for a given gain value and number of
% multiplication pixels.
mean_gain = 100;
mult_prob = mean_gain.^(1 / num_mult) - 1;

% Instantiate parallel and serial sections of the register
test_det.init_par_sec(fwc_par, cte, dcr, cicr_par, wl_ref_qe, iqe_ref);
test_det.init_srl_reg(fwc_ser, cte, dcr, cicr_srl, mult_prob, wl_ref_qe, ...
    iqe_ref);

% Instantiate the controller
cam_gain = 21.9739974;
adc_bits = 14;
adc_offset = 0;
vert_freq = 0.9e6;
horz_freq = 1e6;

test_cntrl = Controller(cam_gain, adc_bits, adc_offset, test_det, ...
    vert_freq, horz_freq);

% Damage the detector
% test_cntrl.place_hot_dark_par(0.05, dcr * 10);
% test_cntrl.place_cti_par(0.05, 0.9);
% test_cntrl.place_hot_dark_srl(0.05, dcr * 10);
% test_cntrl.place_cti_srl(0.05, 0.9);

% Set up the observation parameters
t_exp = 100; % s - exposure time
wl_ob = [1 10 215 300 600 1000];

% Set up random soft x-ray background
ph_map_1nm = 0.001 * rand(par_len, par_wid);

ph_map_10nm = 1 * [...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 1 0 0 0 0 0 1 0 0 0 0 0; ...
    0 0 1 0 0 0 1 0 0 0 0 0 0; ...
    0 0 0 1 0 1 0 0 0 0 0 0 0; ...
    0 0 0 0 1 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0;];

ph_map_215nm = 1 * [...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 1 1 1 1 0; ...
    0 0 0 0 0 0 0 0 0 1 0 0 0; ...
    0 0 0 0 0 0 0 0 0 1 0 0 0; ...
    0 0 0 0 0 0 0 0 0 1 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0;];

ph_map_300nm = 1 * [...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 1 0 0 0 1 0 0 0 0 0 0 0; ...
    0 1 1 0 0 1 0 0 0 0 0 0 0; ...
    0 1 0 1 0 1 0 0 0 0 0 0 0; ...
    0 1 0 0 1 1 0 0 0 0 0 0 0; ...
    0 1 0 0 0 1 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0;];

ph_map_600nm = 1 * [...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 1 1 1 0 0 0; ...
    0 0 0 0 0 0 0 1 0 0 0 0 0; ...
    0 0 0 0 0 0 0 1 1 1 0 0 0; ...
    0 0 0 0 0 0 0 0 0 1 0 0 0; ...
    0 0 0 0 0 0 0 1 1 1 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0;];

ph_map_1000nm = 1 * [...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 0 0 0 0 0;];

ph_map_comp = cat(3, ph_map_1nm, ph_map_10nm, ph_map_215nm, ...
    ph_map_300nm, ph_map_600nm, ph_map_1000nm);

% Plot the input photon maps
figure();
subplot(3, 2, 1);
imagesc(ph_map_1nm);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'photons s^{-1}';
ylabel('Vertical Position, Pixels');
title('1 nm Input');

subplot(3, 2, 2);
imagesc(ph_map_10nm);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'photons s^{-1}';
title('10 nm Input');

subplot(3, 2, 3);
imagesc(ph_map_215nm);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'photons s^{-1}';
ylabel('Vertical Position, Pixels');
title('215 nm Input');

subplot(3, 2, 4);
imagesc(ph_map_300nm);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'photons s^{-1}';
title('300 nm Input');

subplot(3, 2, 5);
imagesc(ph_map_600nm);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'photons s^{-1}';
xlabel('Horizontal Position, Pixels');
ylabel('Vertical Position, Pixels');
title('600 nm Input');

subplot(3, 2, 6);
imagesc(ph_map_1000nm);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'photons s^{-1}';
xlabel('Horizontal Position, Pixels');
title('1000 nm Input');

% Perform the exposure
if read_opt == 1
    test_cntrl.sng_int_full_fr_em(ph_map_comp, wl_ob, t_exp, 1);
else
    test_cntrl.sng_int_full_fr_std(ph_map_comp, wl_ob, t_exp, 1);
end

rtrn_img = test_cntrl.apply_cam_gain();

% Show the output
% Output in DN
figure();
imagesc(rtrn_img);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'DN';
set(gca, 'ColorScale', 'log');
xlabel('Horizontal Position, Pixels');
ylabel('Vertical Position, Pixels');
title('Final Image in DN');

% Output in e-
figure();
imagesc(test_cntrl.image_store);
axis image;
c = colorbar('location', 'eastoutside');
c.Label.String = 'e^-';
set(gca, 'ColorScale', 'log');
xlabel('Horizontal Position, Pixels');
ylabel('Vertical Position, Pixels');
title('Final Image in e^-');