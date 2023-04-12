%% Nicholas Jones - njones31@vt.edu
% Driver script for the detector. Characterizes CIC for a simulated EMCCD.
% Must have estimates for the read noise.
close all;
clear;
clc;

%% Control loading of saved data
load_data = 0;
save_folder = 'char_data_archive\';

if load_data
    data_file = 'sim_cic_data_20230310225252.mat';
    load([save_folder data_file]);
else
    parpool(4);
end

%% If not loading data, setup the EMCCD simulation and collect bias data
if ~load_data
    % Utilize an example iQE curve
    wl_ref_qe = 100 : 1 : 1100; % nm - wavelength vector
    iqe_ref = 0.9 * ones(size(wl_ref_qe)); % Example iQE (0 - 1)
    
    % Define detector properties
    par_len = 25;               % Pixels - Parallel section length
    par_wid = 25;               % Pixels - Parallel section width
    num_std_os = 16;            % Pixels - overscan elelemtns
    num_ad_elem = 16;           % Pixels - standard and corner elements
    num_mult = 50;              % Pixels - multiplication elements
    em_read_noise = 92.798;     % e- - read noise rms on the em register
    em_read_bias = 24030.448;   % e- - bias on the em register
    std_read_noise = 10;        % e- - read noise rms on the standard
                                % register
    std_read_bias = 15;         % e- - bias on the standard register
    
    % Define pixel properties
    fwc_par = 80000;    % e- - Parallel pixel well depth
    fwc_ser = 730000;   % e- - Horizontal pixel well depth
    cte = 0.999989;     % Charge Transfer Efficiency     
    dcr = 7.0685e-4;    % e- pixel^-1 s^-1 - Dark current rate
    
    cic_rate = 0.0031;          % e- pixel^-1 frame^-1 - CIC rate.
    cicr_par = cic_rate / 92;   % e- pixel^-1 transfer^-1 - CIC rate. NOT
                                % per frame.
    cicr_srl = cic_rate / 92;   % e- pixel^-1 transfer^-1 - CIC rate. NOT
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
    % Get inputs for threshold and read noise values for photon-counting
    cic_thresh = input('Enter the photon-counting threshold multiplier: ');
    cic_read_est = input('Enter the read noise estimate in DN: ');
    cic_cam_gain_est = input('Enter the camera gain estimate in DN: ');

    t_exp_c = 0;                                % Bias exposure time
    wl_ob_c = 600;                              % nm - Bias exposure
                                                % wavelength
    ph_map_600nm_c = zeros(par_len, par_wid);   % photons - Bias exposure
                                                % photon map. Set to zero.
    num_cic = 1000;                             % Number of bias frames to
                                                % collect
    
    cic_cube = zeros(par_len, par_wid, num_cic);
    
    % Perform the exposures
    parfor i = 1 : num_cic
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
        disp(['Performing CIC exposure: ' num2str(i)]);

        test_cntrl.sng_int_full_fr_em(ph_map_600nm_c, wl_ob_c, t_exp_c, ...
            1, 1);
        cic_cube(:, :, i) = test_cntrl.apply_cam_gain();

        disp(['Time to complete ' num2str(i) ': ' num2str(toc) ' s']);
    end
end

% Find the mean bias frame
cic_bias = median(cic_cube, 3);
cic_std = std(cic_cube, 1, 3);

% Remove bias from the CIC frames
cic_cube_br = cic_cube - cic_bias;

% Display a histogram of the data
figure();
histogram(cic_cube_br, 'BinWidth', 1);
set(gca, 'YScale', 'Log');
xlim([-Inf 25 * cic_read_est]);
xline(6 * cic_read_est, 'r:', 'Label', '6\sigma_R');
xline(20 * cic_read_est, 'r:', 'Label', '20\sigma_R');
xlabel('DN');
ylabel('Counts');
title('Simulated CIC Data Histogram (Bias Subtracted)');

figure();
imagesc(cic_bias);
colormap gray;
axis image;
xlabel('Pixels');
ylabel('Pixels');
c = colorbar();
c.Label.String = 'DN';
title('Simulated Median CIC Bias Frame');

figure();
imagesc(cic_std);
axis image;
xlabel('Pixels');
ylabel('Pixels');
c = colorbar();
c.Label.String = 'DN';
title('Simulated Standard Deviation of CIC Bias Frame');

% Attempt to determine the CIC using photon counting
cic_cube_pc = cic_cube_br > (cic_thresh * cic_read_est);

cic_pc_frame_count = squeeze(sum(cic_cube_pc, [1 2]));

% Determine the CIC rate in e- pix^-1 fr^-1
cic_rate_est_pc = mean(cic_pc_frame_count) ./ (par_len * par_wid);

% Fit an exponential probability density function to the data between 6 and
% 20 times the read noise. Will determine the mean gain as well.
fit_data = cic_cube_br(cic_cube_br >= 6 * cic_read_est & cic_cube_br ...
    <= 20 * cic_read_est);

[gain_fit, g_ci] = expfit(fit_data);

gain_fit_est = gain_fit * cic_cam_gain_est;

cic_count_est = length(fit_data) / (exp(-20 * cic_read_est / gain_fit) *...
    (-1 + exp(14 * cic_read_est / gain_fit)));

% Determine the CIC rate in e- pix^-1 fr^-1
cic_rate_est_fit = cic_count_est ./ (num_cic * par_len * par_wid);

%% Save the data to a .mat file
if ~load_data
    file_name = [save_folder 'sim_cic_data_' datestr(datetime('now'), ...
            'yyyymmddHHMMSS')];
    
    save(file_name, 'wl_ref_qe', 'iqe_ref', 'par_len', 'par_wid', ...
        'num_std_os', 'num_ad_elem', 'num_mult', 'em_read_noise', ...
        'em_read_bias', 'std_read_noise', 'std_read_bias', 'fwc_par', ...
        'fwc_ser', 'cte', 'dcr', 'cic_rate', 'cicr_par', 'cicr_srl', ...
        'mean_gain', 'cam_gain', 'adc_bits', 'adc_offset', 'vert_freq', ...
        'horz_freq', 'cic_thresh', 'cic_read_est', 'cic_cam_gain_est', ...
        't_exp_c', 'wl_ob_c','ph_map_600nm_c', 'num_cic', 'cic_cube', ...
        '-v7.3');

    delete(gcp('nocreate'));
end