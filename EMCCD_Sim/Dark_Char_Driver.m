%% Nicholas Jones - njones31@vt.edu
% Driver script for the detector. Characterizes dark current for a
% simulated EMCCD. Must have estimates for the read noise and CIC.
close all;
clear;
clc;

%% Control loading of saved data
load_data = 1;
save_folder = 'char_data_archive\';

if load_data
    data_file = 'sim_dark_data_20230311155334.mat';
    load([save_folder data_file]);
else
    parpool(8);
end

%% If not loading data, setup the EMCCD simulation and collect bias data
if ~load_data
    % Utilize an example iQE curve
    wl_ref_qe = 100 : 1 : 1100; % nm - wavelength vector
    iqe_ref = 0.9 * ones(size(wl_ref_qe)); % Example iQE (0 - 1)
    
    % Define detector properties
    par_len = 50;               % Pixels - Parallel section length
    par_wid = 50;               % Pixels - Parallel section width
    num_std_os = 16;            % Pixels - overscan elelemtns
    num_ad_elem = 16;           % Pixels - standard and corner elements
    num_mult = 50;              % Pixels - multiplication elements
    em_read_noise = 92.798;     % e- - read noise rms on the em register
    em_read_bias = 23657.194;   % e- - bias on the em register
    std_read_noise = 10;        % e- - read noise rms on the standard
                                % register
    std_read_bias = 15;         % e- - bias on the standard register
    
    % Define pixel properties
    fwc_par = 80000;    % e- - Parallel pixel well depth
    fwc_ser = 730000;   % e- - Horizontal pixel well depth
    cte = 0.999989;     % Charge Transfer Efficiency     
    dcr = 7.0685e-4;    % e- pixel^-1 s^-1 - Dark current rate
    
    cic_rate = 0.0031;          % e- pixel^-1 frame^-1 - CIC rate.
    cicr_par = cic_rate / 117;  % e- pixel^-1 transfer^-1 - CIC rate. NOT
                                % per frame.
    cicr_srl = cic_rate / 117;  % e- pixel^-1 transfer^-1 - CIC rate. NOT
                                % per frame.
    
    % Calculate multiplication probability for a given gain value and
    % number of multiplication pixels.
    mean_gain = 1000;
    
    % Define controller properties
    cam_gain = 17.7740;
    adc_bits = 16;
    adc_offset = 0;
    vert_freq = 1e6;
    horz_freq = 10e6;
end

% Plot the example iQE curve
figure();
plot(wl_ref_qe, iqe_ref);
xlim([wl_ref_qe(1) wl_ref_qe(end)]);
ylim([0 1]);
xlabel('Wavelegnth, nm');
ylabel('iQE');

%% CIC Data Collection
% If not loading data, simulate bias data collection
if ~load_data
    dark_thresh = input(...
        'Enter the photon-counting threshold multiplier: ');
    dark_read_est = input('Enter the read noise estimate in DN: ');
    dark_cam_gain_est = input('Enter the camera gain estimate in DN: ');
    dark_cic_est = input('Enter the CIC estimate in e-: ');

    t_exp_d = 100;                              % s - dark exposure time
    wl_ob_d = 600;                              % nm - dark exposure
                                                % wavelength
    ph_map_600nm_d = zeros(par_len, par_wid);   % photons - dark exposure
                                                % photon map. Set to zero.
    num_dark = 144;                             % Number of dark frames to
                                                % collect
    
    dark_cube = zeros(par_len, par_wid, num_dark);
    
    % Perform the exposures
    parfor i = 1 : num_dark
        % Instantiate the detector
        test_det = Detector(par_wid, par_len, num_std_os, num_ad_elem, ...
            num_mult, em_read_noise, em_read_bias, std_read_noise, ...
            std_read_bias);

        % Instantiate parallel and serial sections of the register
        test_det.init_par_sec(fwc_par, cte, dcr, cicr_par, wl_ref_qe, ...
            iqe_ref);
        test_det.init_srl_reg(fwc_ser, cte, dcr, cicr_srl, mean_gain, ...
            wl_ref_qe, iqe_ref);

        % Instantiate the controller
        test_cntrl = Controller(cam_gain, adc_bits, adc_offset, ...
            test_det, vert_freq, horz_freq);

        tic;
        disp(['Performing Dark exposure: ' num2str(i)]);

        test_cntrl.sng_int_full_fr_em(ph_map_600nm_d, wl_ob_d, t_exp_d, ...
            1, 1);
        dark_cube(:, :, i) = test_cntrl.apply_cam_gain();

        disp(['Time to complete ' num2str(i) ': ' num2str(toc) ' s']);
    end
end

% Find the mean bias frame
dark_bias = median(dark_cube, 3);
dark_std = std(dark_cube, 1, 3);

% Remove bias from the CIC frames
dark_cube_br = dark_cube - dark_bias;

% Display a histogram of the data
figure();
histogram(dark_cube_br, 'BinWidth', 1);
set(gca, 'YScale', 'Log');
xlabel('DN');
ylabel('Counts');
title('Simulated Dark Data Histogram (Bias Subtracted)');

figure();
imagesc(dark_bias);
colormap gray;
axis image;
xlabel('Pixels');
ylabel('Pixels');
c = colorbar();
c.Label.String = 'DN';
title('Simulated Median Dark Bias Frame');

figure();
imagesc(dark_std);
axis image;
xlabel('Pixels');
ylabel('Pixels');
c = colorbar();
c.Label.String = 'DN';
title('Simulated Standard Deviation of Median Dark Bias Frame');

% Attempt to determine the dark current using photon counting
dark_cube_pc = dark_cube_br > (dark_thresh * dark_read_est);

dark_pc_frame_count = squeeze(sum(dark_cube_pc, [1 2]));
dark_rate_est_pc = (mean(dark_pc_frame_count) ./ (par_len * par_wid) - ...
    dark_cic_est) / t_exp_d; % e- pix^-1 s^-1

% Fit an exponential probability density function to the data between 5 and
% 20 times the read noise. Will determine the mean gain as well. Assumes
% the data is well descibed by having < 1 input electron at the beginning
% of charge multiplication. Based on the procedure in 'Extreme Faint Flux
% Imaging with an EMCCD' by Olivier Daigle et al. 2009.
fit_data = dark_cube_br(dark_cube_br >= 6 * dark_read_est & dark_cube_br...
    <= 20 * dark_read_est);

[gain_fit, g_ci] = expfit(fit_data);

gain_fit_est = gain_fit * dark_cam_gain_est;

dark_count_est = length(fit_data) / (exp(-20 * dark_read_est / ...
    gain_fit) * (-1 + exp(14 * dark_read_est / gain_fit)));

% Determine the CIC rate in e- pix^-1 fr^-1
dark_rate_est_fit = ((dark_count_est ./ (num_dark * par_len * par_wid))...
    - dark_cic_est) / t_exp_d;

%% Save the data to a .mat file
if ~load_data
    file_name = [save_folder 'sim_dark_data_' datestr(datetime('now'), ...
            'yyyymmddHHMMSS') '.mat'];
    
    save(file_name, 'wl_ref_qe', 'iqe_ref', 'par_len', 'par_wid', ...
        'num_std_os', 'num_ad_elem', 'num_mult', 'em_read_noise', ...
        'em_read_bias', 'std_read_noise', 'std_read_bias', 'fwc_par', ...
        'fwc_ser', 'cte', 'dcr', 'cic_rate', 'cicr_par', 'cicr_srl', ...
        'mean_gain', 'cam_gain', 'adc_bits', 'adc_offset', 'vert_freq', ...
        'horz_freq', 'dark_thresh', 'dark_read_est', ...
        'dark_cam_gain_est', 'dark_cic_est', 't_exp_d', 'wl_ob_d', ...
        'ph_map_600nm_d', 'num_dark', 'dark_cube', '-v7.3');

    delete(gcp('nocreate'));
end