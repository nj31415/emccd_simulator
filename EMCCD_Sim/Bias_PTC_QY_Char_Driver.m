%% Nicholas Jones - njones31@vt.edu
% Driver script for the detector. Performs a photon transfer curve
% operation, varying the integration time for the exposure to produce
% different intensity measurements. Also produces bias frames.
close all;
clear;
clc;

%% Control loading of saved data
load_data = 0;
save_folder = 'char_data_archive\';

if load_data
    data_file = 'sim_bias_ptc_600_data_20230309135046.mat';
    load([save_folder data_file]);
else
    parpool(8);
end

%% If not loading data, setup the EMCCD simulation and collect bias data
if ~load_data
    % Utilize an example iQE curve
    wl_ref_qe = 0.21 : 1 : 1100; % nm - wavelength vector
    iqe_ref = 0.9 * ones(size(wl_ref_qe)); % Example iQE (0 - 1)
    
    % Define detector properties
    par_len = 50;               % Pixels - Parallel section length
    par_wid = 50;               % Pixels - Parallel section width
    num_std_os = 16;            % Pixels - overscan elelemtns
    num_ad_elem = 16;           % Pixels - standard and corner elements
    num_mult = 50;              % Pixels - multiplication elements
    em_read_noise = 92.798;     % e- - read noise rms on the em register
    em_read_bias = 33397.346; % e- - bias on the em register
    std_read_noise = 10;        % e- - read noise rms on the standard
                                % register
    std_read_bias = 15;         % e- - bias on the standard register
    
    % Define pixel properties
    fwc_par = 80000;    % e- - Parallel pixel well depth
    fwc_ser = 730000;   % e- - Horizontal pixel well depth
    cte = 0.999989;     % Charge Transfer Efficiency     
    dcr = 7.0685e-4;    % e- pixel^-1 s^-1 - Dark current rate
    
    cic_rate = 0.0031;          % e- pixel^-1 frame^-1 - CIC rate.

    % To translate the CIC rate from per frame to per transfer, the CIC
    % rate is divided by the average number of transfers. The script
    % cic_fr_2_tr.m can help determine this number.
    cicr_par = cic_rate / 117;  % e- pixel^-1 transfer^-1 - CIC rate. NOT
                                % per frame.
    cicr_srl = cic_rate / 117;  % e- pixel^-1 transfer^-1 - CIC rate. NOT
                                % per frame.
    
    % Calculate multiplication probability for a given gain value and
    % number of multiplication pixels.
    mean_gain = 1;
    
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

%% Bias Data Collection
% If not loading data, simulate bias data collection
if ~load_data
    t_exp_b = 0;                                % Bias exposure time
    wl_ob_b = 600;                              % nm - Bias exposure
                                                % wavelength
    ph_map_600nm_b = zeros(par_len, par_wid);   % photons - Bias exposure
                                                % photon map. Set to zero.
    num_bias = 50;                             % Number of bias frames to
                                                % collect
    
    bias_cube = zeros(par_len, par_wid, num_bias);
    
    % Perform the exposures
    parfor i = 1 : num_bias
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
        test_cntrl = Controller(cam_gain, adc_bits, adc_offset, test_det, ...
        vert_freq, horz_freq);

        disp(['Performing bias exposure: ' num2str(i)]);

        test_cntrl.sng_int_full_fr_em(ph_map_600nm_b, wl_ob_b, t_exp_b, ...
            1, 1);
        bias_cube(:, :, i) = test_cntrl.apply_cam_gain();
    end
end

% Display a histogram of the data
figure();
histogram(bias_cube, 'BinWidth', 1);
set(gca, 'YScale', 'Log');
xlabel('DN');
ylabel('Counts');
title('Simulated Bias Data Histogram');

% Find the mean bias frame
bias_median = median(bias_cube, 3);
bias_std = std(bias_cube, 1, 3);

figure();
imagesc(bias_median);
colormap gray;
axis image;
xlabel('Pixels');
ylabel('Pixels');
c = colorbar();
c.Label.String = 'DN';
title('Simulated Median Bias Frame');

figure();
imagesc(bias_std);
axis image;
xlabel('Pixels');
ylabel('Pixels');
c = colorbar();
c.Label.String = 'DN';
title('Simulated Standard Deviation of Bias Frame');

%% Photon Transfer Curve
% If not loading data, simulate PTC data generation
if ~load_data
    % Set up the PTC observation parameters
    num_ptc = 40;
    t_exp_ptc = logspace(1, 5, num_ptc); % s - exposure time
    wl_ob_ptc = 0.21; % nm - light wavelength in nm
    
    ph_map_ptc = ones(par_len, par_wid);
    
    % Set up output variables
    t_vec = zeros(size(t_exp_ptc)); % Track how long each iteration takes
    
    ptc_cube = zeros(par_len, par_wid, num_ptc, 2);
    
    % Perform the exposures
    parfor i = 1 : num_ptc
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
        test_cntrl = Controller(cam_gain, adc_bits, adc_offset, test_det, ...
        vert_freq, horz_freq);

        tic;
        disp(['Performing PTC exposure: ' num2str(i)]);

        % First exposure. Blanks array before integration
        test_cntrl.sng_int_full_fr_em(ph_map_ptc, wl_ob_ptc, ...
            t_exp_ptc(i), 1, 1);
        rtrn_img_1 = test_cntrl.apply_cam_gain();
        
        % Second exposure. Does not blank array before integration
        test_cntrl.sng_int_full_fr_em(ph_map_ptc, wl_ob_ptc, ...
            t_exp_ptc(i), 0, 1);
        rtnr_img_2 = test_cntrl.apply_cam_gain();

        ptc_cube(:, :, i, :) = cat(4, rtrn_img_1, rtnr_img_2);
        
        t_vec(i) = toc;
    end
end

ptc_mean_vec = squeeze(mean(ptc_cube - bias_median, [1, 2, 4]));
ptc_diff_noise = squeeze(std(ptc_cube(:, :, :, 1) - ...
    ptc_cube(:, :, :, 2), 1, [1 2]) ./ sqrt(2));

% Sort PTC vectors in terms of increasing signal level
[ptc_mean_vec, sort_idx] = sort(ptc_mean_vec);
ptc_diff_noise = ptc_diff_noise(sort_idx);

% Remove negative signal levels
n_idx = ptc_mean_vec < 0;
sort_idx(n_idx) = [];
ptc_mean_vec(n_idx) = [];
ptc_diff_noise(n_idx) = [];

figure();
loglog(ptc_mean_vec, ptc_diff_noise, 'k*');
xlabel('Signal, S(DN)');
ylabel('Noise, \sigma_S(DN)');

figure();
semilogx(t_exp_ptc, t_vec / 2, 'k*');
xlabel('Exposure Time, s');
ylabel('Computation Time, s');

% Display the PTC vectors
disp('      Signal | Noise');
disp([ptc_mean_vec ptc_diff_noise]);

drawnow;

% Estimate the read noise
num_rn_regime = input('Number of data points in the read noise regime: ');
sigma_rn_est = mean(ptc_diff_noise(1 : num_rn_regime));
sigma_rn_std = std(ptc_diff_noise(1 : num_rn_regime));
disp(['Read noise estimate (DN): ' num2str(sigma_rn_est)]);
disp(['in e-: ' num2str(sigma_rn_est * cam_gain)]);
disp(['Standard deviation (DN): ' num2str(sigma_rn_std)]);
disp([' in e-: ' num2str(sigma_rn_std * cam_gain)]);

% Attempt to calculate the camera gain constant
num_sat_elem = input('Number of saturated data points: ');
num_ex_elem = input('Number of data points to exclude at beginning: ');
lin_idx = num_ex_elem + 1 : length(ptc_mean_vec) - num_sat_elem;
gain_cnst_est_vec = ptc_mean_vec(lin_idx) ./ ...
    (ptc_diff_noise(lin_idx).^2 - sigma_rn_est.^2);
gain_cnst_est = mean(gain_cnst_est_vec);
gain_cnst_std = std(gain_cnst_est_vec, 1);
disp(['Camera gain constant estimate (e- DN^-1): ' ...
    num2str(gain_cnst_est)]);
disp(['Standard deviation: ' num2str(gain_cnst_std)]);

%% Save the data to a .mat file
if ~load_data
    file_name = [save_folder 'sim_bias_ptc_' num2str(wl_ob_ptc) '_data_'...
        datestr(datetime('now'), 'yyyymmddHHMMSS') '.mat'];
    
    save(file_name, 'wl_ref_qe', 'iqe_ref', 'par_len', 'par_wid', ...
        'num_std_os', 'num_ad_elem', 'num_mult', 'em_read_noise', ...
        'em_read_bias', 'std_read_noise', 'std_read_bias', 'fwc_par', ...
        'fwc_ser', 'cte', 'dcr', 'cic_rate', 'cicr_par', 'cicr_srl', ...
        'mean_gain', 'cam_gain', 'adc_bits', 'adc_offset', 'vert_freq', ...
        'horz_freq', 't_exp_b', 'wl_ob_b', 'ph_map_600nm_b', 'num_bias',...
        'bias_cube', 'num_ptc', 't_exp_ptc', 'wl_ob_ptc', ...
        'ph_map_ptc', 't_vec', 'ptc_cube', '-v7.3');

    delete(gcp('nocreate'));
end