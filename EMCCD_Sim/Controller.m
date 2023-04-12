% Nicholas Jones - njones31@vt.edu
% The controller class is used to coordinate the integration and readout of
% the Detector. The controller class also applies controller effects such
% as the conversion gain 
classdef Controller < handle
    properties (SetAccess = public)
        cam_gain    % e- DN^-1 - camera gain constant
        adc_bits    % # bits - number of bits in the ADC
        adc_offset  % # bits - ADC offset corresponding to 0 voltage / e-
        vert_freq
        horz_freq

        det    % Detector object for the controller
    end

    properties (SetAccess = private)
        image_store % Array for storing the output image
    end
    
    methods
        %% Constructor for a Controller object.
        % cam_gain:     Float, camera gain constant in e- DN^-1
        % adc_bits:     Int, number of bits in the ADC
        % adc_offset:   Int, offset of the zero-level of the ADC
        % detector:     Detector, completely specified detector object to
        %               use with the controller
        % vert_freq:    Float, vertical transfer frequency in Hz
        % horz_freq:    Float, horizontal trasnfer frequency in Hz
        function obj = Controller(cam_gain, adc_bits, adc_offset, ...
                detector, vert_freq, horz_freq)
            obj.cam_gain = cam_gain;
            obj.adc_bits = adc_bits;
            obj.adc_offset = adc_offset;
            obj.det = detector;
            obj.vert_freq = vert_freq;
            obj.horz_freq = horz_freq;
            
            % Set up the image storage array based on detector object
            obj.init_image_store();
        end

        %% Function to get the camera gain constant.
        % Outputs:
        % cam_gain: Float, the camera gain constant in e- DN^-1
        function cam_gain = get.cam_gain(obj)
            cam_gain = obj.cam_gain;
        end

        %% Function set the camera gain constant.
        % Inputs:
        % cam_gain: Float, the camera gain constant in e- DN^-1
        function set.cam_gain(obj, cam_gain)
            obj.cam_gain = cam_gain;
        end
        
        %%  Function to get the number of bits for the ADC.
        % Outputs:
        % adc_bits: Int, number of bits in the ADC
        function adc_bits = get.adc_bits(obj)
            adc_bits = obj.adc_bits;
        end

        %% Function to set the number of bits for the ADC.
        % Inputs:
        % adc_bits: Int, number of bits in the ADC
        function set.adc_bits(obj, adc_bits)
            obj.adc_bits = adc_bits;
        end
        
        %% Function to get the ADC offset
        % Outputs:
        % adc_offset:   Int, ADC offset
        function adc_offset = get.adc_offset(obj)
            adc_offset = obj.adc_offset;
        end

        %% Function to set the ADC offset. Enforces the adc_offset being a 
        % whole number by using the round function.
        % Inputs:
        % adc_offset:   Int, ADC offset
        function set.adc_offset(obj, adc_offset)
            obj.adc_offset = round(adc_offset);
        end
        
        %% Function to get the detector of the controller
        % Outputs:
        % detector: Detector object
        function detector = get.det(obj)
            detector = obj.det;
        end

        %% Function to set the detector of the controller
        % Inputs:
        % detector: Detector object, must be set up prior to being added to
        %           the detector
        function set.det(obj, detector)
            if detector.check_init()
                obj.det = detector;
            else
                error(['Controller::set.detector: Detector has not ' ...
                    'been initiated']);
            end
        end

        %% Function to get the vertical transfer frequency
        % Outputs:
        % vert_freq:   Float, Vertical transfer frequency in Hz
        function vert_freq = get.vert_freq(obj)
            vert_freq = obj.vert_freq;
        end

        %% Function to set the vertical transfer frequency
        % Inputs:
        % vert_freq:   Float, Vertical transfer frequency in Hz
        function set.vert_freq(obj, vert_freq)
            obj.vert_freq = vert_freq;
        end

        %% Function to get the horizontal transfer frequency
        % Outputs:
        % horz_freq:   Float, Horizontal transfer frequency in Hz
        function horz_freq = get.horz_freq(obj)
            horz_freq = obj.horz_freq;
        end

        %% Function to set the horizontal transfer frequency
        % Inputs:
        % horz_freq:   Float, Horizontal transfer frequency in Hz
        function set.horz_freq(obj, horz_freq)
            obj.horz_freq = horz_freq;
        end

        %% Function to get the image from the image storage
        % Outputs:
        % image_store:  Float array, the image data
        function image_store = get.image_store(obj)
            image_store = obj.image_store;
        end

        %% Function to directly set the image store. This method should not
        % be used typically.
        % Inputs:
        % image_store: Float array, image data
        function set.image_store(obj, image_store)
            obj.image_store = image_store;
        end

        %% Function to set up the image storage array based on detector
        % object. The detector for the object must be defined before this
        % function is called
        function init_image_store(obj)
            obj.image_store = zeros(obj.det.par_len, ...
                obj.det.par_wid);
        end

        %% Function to insert data into the image store
        % Inputs:
        % pix_data: Float, data from the detector to place in the image
        %           store in e-
        % x:        Int, width coordinate in the image
        % y:        Int, length coordinate in the image
        function insert_img_data(obj, pix_data, x, y)
            if x < 1 || x > obj.det.par_wid
                error(['Controller::insert_img_data: Invalid x ' ...
                    'coordinate']);
            end

            if y < 1 || y > obj.det.par_len
                error(['Controller::insert_img_data: Invalid y ' ...
                    'coordinate']);
            end

            obj.image_store(y, x) = pix_data;
        end

        %% Function to convert the image from e- to DN
        % Outputs:
        % dn_image: Int array, image_store converted into DN. image_store
        %           is not changed by this function
        function dn_image = apply_cam_gain(obj)
            dn_image = obj.e_2_dn(obj.image_store);
        end

        %% Function to perform a single integration, full frame readout of
        % the detector using the EM register. Image data is stored in the
        % image store of the controller.
        % Inputs:
        % ph_map:   Float array, photon map corresponding to the pixel of
        %           the parallel section of the detector. Photon levels
        %           correspond to the mean photon rate in photons sec^-1.
        %           The third dimension corresponds to different
        %           wavelengths of light
        % wl_vec:   Float vector, list of wavelengths in nm corresponding
        %           to the 3rd dimension of the photon map
        % exp_int:  Float, exposure integration time in sec
        % opt_1:    Boolean, set to 1 to clear the detector before the
        %           exposure
        % opt_2:    Boolean, set to 1 to close the 'shutter' of the camera
        %           during read out. Sets the input photon map to 0 for the
        %           read-out loop Set to zero to keep the shutter open.
        function sng_int_full_fr_em(obj, ph_map, wl_vec, exp_int, opt_1, ...
                opt_2)
            if opt_1
                obj.det.clear_detector();
            end
            
            % Perform an exposure of the detector
            obj.det.integrate(ph_map, wl_vec, exp_int);
            
            % Dump the serial register to prepare for read out
            obj.det.dump_drain(ph_map, wl_vec, 1, obj.det.get_srl_len(),...
                obj.vert_freq, 0);

            % Set up storage variables for detector sizing parameters
            b = obj.det.num_ad_elem;
            c = obj.det.num_mult;
            par_size = obj.det.par_len * obj.det.par_wid;

            % Set up the image index to track where in the read-out the
            % process is.
            img_idx = -(b + c);
            
            % If closing the shutter during readout, blanks the incoming
            % photon map
            if opt_2
                ph_map = zeros(size(ph_map));
            end

            % Enter the read-out loop
            while img_idx < par_size
                if mod(img_idx + b + c, obj.det.par_wid) == 0 && ...
                        img_idx + b + c < par_size
                    obj.det.fw_clock_par(ph_map, wl_vec, obj.vert_freq, 1);
                end

                img_idx = img_idx + 1;

                if img_idx > 0
                    x_coord = obj.det.par_wid - ...
                        mod(img_idx - 1, obj.det.par_wid);
                    y_coord = ceil(img_idx / obj.det.par_wid);

                    obj.insert_img_data(...
                        obj.det.em_read(ph_map, wl_vec, obj.horz_freq), ...
                        x_coord, y_coord);
                else
                    obj.det.em_read(ph_map, wl_vec, obj.horz_freq);
                end
            end
        end

        %% Function to perform a single integration, full frame readout of
        % the detector using the standard register. Image data is stored in
        % the image store of the controller.
        % Inputs:
        % ph_map:   Float array, photon map corresponding to the pixel of
        %           the parallel section of the detector. Photon levels
        %           correspond to the mean photon rate in photons sec^-1.
        %           The third dimension corresponds to different
        %           wavelengths of light
        % wl_vec:   Float vector, list of wavelengths in nm corresponding
        %           to the 3rd dimension of the photon map
        % exp_int:  Float, exposure integration time in sec
        % opt_1:    Boolean, set to 1 to clear the detector before the
        %           exposure
        % opt_2:    Boolean, set to 1 to close the 'shutter' of the camera
        %           during read out. Sets the input photon map to 0 for the
        %           read-out loop Set to zero to keep the shutter open.
        function sng_int_full_fr_std(obj, ph_map, wl_vec, exp_int, opt_1, opt_2)
            if opt_1
                obj.det.clear_detector();
            end
            
            % Perform an exposure of the detector
            obj.det.integrate(ph_map, wl_vec, exp_int);
            
            % Dump the serial register to prepare for read out
            obj.det.dump_drain(ph_map, wl_vec, 1, obj.det.get_srl_len(),...
                obj.vert_freq, 0);

            % Set up storage variables for detector sizing parameters
            a = obj.det.num_std_os;
            par_size = obj.det.par_len * obj.det.par_wid;

            % Set up the image index to track where in the read-out the
            % process is.
            img_idx = -a;

            % If closing the shutter during readout, blanks the incoming
            % photon map
            if opt_2
                ph_map = zeros(size(ph_map));
            end

            % Enter the read-out loop
            while img_idx < par_size
                if mod(img_idx + a, obj.det.par_wid) == 0 && ...
                        img_idx + a < par_size
                    obj.det.fw_clock_par(ph_map, wl_vec, obj.vert_freq, 1);
                end

                img_idx = img_idx + 1;

                if img_idx > 0
                    x_coord = mod(img_idx - 1, obj.det.par_wid) + 1;
                    y_coord = ceil(img_idx / obj.det.par_wid);

                    obj.insert_img_data(...
                        obj.det.std_read(ph_map, wl_vec, obj.horz_freq), ...
                        x_coord, y_coord);
                else
                    obj.det.std_read(ph_map, wl_vec, obj.horz_freq);
                end
            end
        end
        
        %% Function to damage the parallel section of the detector with hot
        % pixels. Damaged pixels are placed randomly.
        % Inputs:
        % frac_pix: Float (0 - 1), fraction of pixels to degrade
        % dcr:      Float, dark current generation rate in e- s^-1
        function place_hot_dark_par(obj, frac_pix, dcr)
            num_pix = floor(obj.det.par_len * obj.det.par_wid * frac_pix);
            rnd_coords = rand(num_pix, 2);
            rnd_coords = [floor(obj.det.par_len * rnd_coords(:, 1)) + 1 ...
                floor(obj.det.par_wid * rnd_coords(:, 2)) + 1];
            
            for i = 1  : num_pix
                obj.det.par_sec(rnd_coords(i, 1), rnd_coords(i, 2)). ...
                    dark_curr_rate = dcr;
            end
        end

        %% Function to randomly degrade CTE in the parallel section.
        % Damaged pixels are placed randomly.
        % Inputs:
        % frac_pix: Float (0 - 1), fraction of pixels to degrade
        % new_cte:  Float (0 - 1), new charge transfer efficiency
        function place_cti_par(obj, frac_pix, new_cte)
            num_pix = floor(obj.det.par_len * obj.det.par_wid * frac_pix);
            rnd_coords = rand(num_pix, 2);
            rnd_coords = [floor(obj.det.par_len * rnd_coords(:, 1)) + 1 ...
                floor(obj.det.par_wid * rnd_coords(:, 2)) + 1];
            
            for i = 1  : num_pix
                obj.det.par_sec(rnd_coords(i, 1), rnd_coords(i, 2)). ...
                    charge_transfer_eff = new_cte;
            end
        end

        %% Function to damage the serial section of the detector with hot
        % pixels. Damaged pixels are placed randomly.
        % Inputs:
        % frac_pix: Float (0 - 1), fraction of pixels to degrade
        % dcr:      Float, dark current generation rate in e- s^-1
        function place_hot_dark_srl(obj, frac_pix, dcr)
            num_pix = floor(obj.det.get_srl_len() * frac_pix);
            rnd_coords = rand(num_pix, 1);
            rnd_coords = floor(obj.det.get_srl_len() * rnd_coords) + 1;
            
            for i = 1  : num_pix
                obj.det.srl_reg(rnd_coords(i)).dark_curr_rate = dcr;
            end
        end

        %% Function to randomly degrade CTE in the serial section. Damaged 
        % pixels are placed randomly.
        % Inputs:
        % frac_pix: Float (0 - 1), fraction of pixels to degrade
        % new_cte:  Float (0 - 1), new charge transfer efficiency
        function place_cti_srl(obj, frac_pix, new_cte)
            num_pix = floor(obj.det.get_srl_len() * frac_pix);
            rnd_coords = rand(num_pix, 1);
            rnd_coords = floor(obj.det.get_srl_len() * rnd_coords) + 1;
            
            for i = 1  : num_pix
                obj.det.srl_reg(rnd_coords(i)).charge_transfer_eff = ...
                    new_cte;
            end
        end

        %% Function to convert from measured pixel value in e- to DN. 
        % Assumes 0 e- corresponds to 0 bits
        % Inputs:
        % e_m:  Float, number of measured e-, in e-
        % Outputs:
        % dn_m: Int, output ADC value, in digital numbers
        function dn_m = e_2_dn(obj, e_m)
            dn_m = round(obj.cam_gain^-1 * e_m) + obj.adc_offset;
            
            % Enforce ADC boundaries (outputs are between 0 and maximum ADC
            % value, given by 2^N_bits - 1).
            dn_m(dn_m > (2^obj.adc_bits - 1)) = (2^obj.adc_bits - 1);
            dn_m(dn_m < 0) = 0;
        end
    end
end

