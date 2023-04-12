%% Nicholas Jones - njones31@vt.edu
% The Detector class is used to represent a detector. The detector
% implements the specific charge trasnfer operations that make up a read
% out by coordinating charge movement through two arrays of pixels. One
% array is the parallel section, where incoming photons are able to
% genearte charge in the array. The other is the serial register, which
% transports the charge to the detector output.

classdef Detector < handle
    properties (SetAccess = private)
        % Define parallel array properties
        par_sec = Pixel.empty();
        par_wid
        par_len
        
        % Define serial register array properties
        srl_reg = Pixel.empty();
        num_std_os
        num_ad_elem
        num_mult
    end

    properties (SetAccess = public)
        em_read_noise
        em_read_bias
        std_read_noise
        std_read_bias
        mean_gain
    end
    
    methods
        %% Constructor for a Detector object. Note this will only set up
        % the basic sizing parameters. Additional functions must be called
        % to instantiate the parallel and serial register pixel arrays with
        % default pixel values.
        % Inputs:
        % par_wid:          Int, width of the parallel section in pixels
        % par_len:          Int, length of the parallel section in pixels
        % num_std_os:       Int, number of overscan elements in the serial
        %                   register before the standard read-out amplifier
        % num_ad_elem:      Int, number of additional elements in the
        %                   serial register before the EM elements
        % num_mult:         Int, number of electron multiplying elements
        % em_read_noise:    Float, rms read-out noise for the EM read-out
        %                   amplifier
        % em_read_bias:     Float, EM read-out bias in e-
        % std_read_noise:   Float, rms read-out noise for the standard 
        %                   read-out amplifier
        % std_read_bias:    Int, standard read-out bias in e-
        function obj = Detector(par_wid, par_len, num_std_os, ...
                num_ad_elem, num_mult, em_read_noise, em_read_bias, ...
                std_read_noise, std_read_bias)
            obj.par_wid = par_wid;
            obj.par_len = par_len;
            obj.num_std_os = num_std_os;
            obj.num_ad_elem = num_ad_elem;
            obj.num_mult = num_mult;
            obj.em_read_noise = em_read_noise;
            obj.em_read_bias = em_read_bias;
            obj.std_read_noise = std_read_noise;
            obj.std_read_bias = std_read_bias;
        end

        %% Function to initialize the parallel section with default pixel
        % properties. Assumes parallel section is non-electron-multiplying.
        % Inputs:
        % full_well_capacity:   e- - Full well capacity
        % charge_transfer_eff:  Charge transfer efficiency (0 - 1)
        % dark_curr_rate:       e- pix^-1 sec^-1 - Dark current
        %                       generation rate
        % cic_rate:             e- pix^-1 tr^-1 - CIC generation rate
        % wl_ref_iqe:           nm - Reference wavelength vector for iQE
        % i_quauntum_eff:       Interacting quantum efficiency (0 - 1).
        %                       Note that this is not the same as quantum
        %                       efficiency
        function init_par_sec(obj, full_well_capacity, ...
                charge_transfer_eff, dark_curr_rate, cic_rate, ...
                wl_ref_iqe, i_quantum_eff)
            % See discussion on initializing handle object arrays at 
            % https://www.mathworks.com/help/matlab/matlab_oop/initializing
            % -arrays-of-handle-objects.html
            obj.par_sec(obj.par_len, obj.par_wid) = ...
                Pixel(full_well_capacity, charge_transfer_eff, ...
                dark_curr_rate, cic_rate, 0, {wl_ref_iqe}, ...
                {i_quantum_eff});

            % Populate pixel properties
            [obj.par_sec.full_well_capacity] = deal(full_well_capacity);
            [obj.par_sec.charge_transfer_eff] = deal(charge_transfer_eff);
            [obj.par_sec.dark_curr_rate] = deal(dark_curr_rate);
            [obj.par_sec.cic_rate] = deal(cic_rate);
            [obj.par_sec.mult_prob] = deal(0);
            [obj.par_sec.i_quantum_eff] = deal({i_quantum_eff});
            [obj.par_sec.wl_ref_iqe] = deal({wl_ref_iqe});
        end

        %% Function to initialize the serial register section with default
        % pixel properties. Assumes constant CIC and dark current on the
        % serial register and that mult_prob is zero except for the
        % electron-multiplying elements. Quantum efficiency can be set as
        % desired by the user, however in normal readout modes the serial
        % register is never exposed to light, so these values will be
        % unused.
        % Inputs:
        % full_well_capacity:   e- - Full well capacity
        % charge_transfer_eff:  Charge transfer efficiency (0 - 1)
        % dark_curr_rate:       e- pix^-1 sec^-1 - Dark current
        %                       generation rate
        % cic_rate:             e- pix^-1 tr^-1 - CIC generation rate
        % mean_gain:            mean gain of the EM register
        % wl_ref_iqe:           nm - Reference wavelength vector for iQE
        % i_quauntum_eff:       Interacting quantum efficiency (0 - 1).
        %                       Note that this is not the same as quantum
        %                       efficiency
        function init_srl_reg(obj, full_well_capacity, ...
                charge_transfer_eff, dark_curr_rate, cic_rate, ...
                mean_gain, wl_ref_iqe, i_quantum_eff)
            % See discussion on initializing handle object arrays at 
            % https://www.mathworks.com/help/matlab/matlab_oop/initializing
            % -arrays-of-handle-objects.html
            obj.srl_reg(obj.get_srl_len()) = Pixel(full_well_capacity, ...
                charge_transfer_eff, dark_curr_rate, cic_rate, 0, ...
                {wl_ref_iqe}, {i_quantum_eff});
            
            % Instantiate the multiplication probability over the entire EM
            % Register
            [obj.srl_reg.mult_prob] = deal(0);

            % Set the mean gain. Will also set the multiplication
            % probability in the EM register to the proper value.
            obj.mean_gain = mean_gain;
            
            % Populate pixel properties
            [obj.srl_reg.full_well_capacity] = deal(full_well_capacity);
            [obj.srl_reg.charge_transfer_eff] = deal(charge_transfer_eff);
            [obj.srl_reg.dark_curr_rate] = deal(dark_curr_rate);
            [obj.srl_reg.cic_rate] = deal(cic_rate);
            [obj.srl_reg.i_quantum_eff] = deal({i_quantum_eff});
            [obj.srl_reg.wl_ref_iqe] = deal({wl_ref_iqe});
        end

        %% Function to get the parallel section width
        % Outputs:
        % par_wid:  Int, parallel section width in pixels
        function par_wid = get.par_wid(obj)
            par_wid = obj.par_wid;
        end

        %% Function to set the parallel section width
        % Inputs:
        % par_wid:  Int, parallel section width in pixels
        function set.par_wid(obj, par_wid)
            obj.par_wid = par_wid;
        end

        %% Function to get the parallel section length
        % Outputs:
        % par_len:  Int, parallel section length in pixels
        function par_len = get.par_len(obj)
            par_len = obj.par_len;
        end

        %% Function to set the parallel section length
        % Inputs:
        % par_len:  Int, parallel section length in pixels
        function set.par_len(obj, par_len)
            obj.par_len = par_len;
        end

        %% Function to get the number of overscan elements before the
        % standard output amplifier
        % Outputs:
        % num_std_os:   Int, number of overscan elements before the
        %               standard output
        function num_std_os = get.num_std_os(obj)
            num_std_os = obj.num_std_os;
        end

        %% Function to set the number of overscan elements before the
        % standard output amplifier
        % Inputs:
        % num_std_os:   Int, number of overscan elements before the
        %               standard output
        function set.num_std_os(obj, num_std_os)
            obj.num_std_os = num_std_os;
        end

        %% Function to get the number of additional elements before the EM
        % register
        % Outputs:
        % num_ad_elem:  Int, number of additional elements before the EM
        %               register
        function num_ad_elem = get.num_ad_elem(obj)
            num_ad_elem = obj.num_ad_elem;
        end

        %% Function to set the number of additional elements before the EM
        % register
        % Inputs:
        % num_ad_elem:  Int, number of additional elements before the EM
        %               register
        function set.num_ad_elem(obj, num_ad_elem)
            obj.num_ad_elem = num_ad_elem;
        end

        %% Function to get the number of multiplication elements
        % Outputs:
        % num_mult: Int, number of multiplication elements
        function num_mult = get.num_mult(obj)
            num_mult = obj.num_mult;
        end

        %% Function to set the number of multiplication elements
        % Inputs:
        % num_mult:    Int, number of multiplication elements
        function set.num_mult(obj, num_mult)
            obj.num_mult = num_mult;
        end

        %% Function to get the rms EM read-out noise
        % Outputs:
        % read_out_noise:   Float, rms EM read-out noise in e-
        function em_read_noise = get.em_read_noise(obj)
            em_read_noise = obj.em_read_noise;
        end

        %% Function to set the rms EM read-out noise
        % Inputs:
        % read_out_noise:   Float, rms EM read-out noise in e-
        function set.em_read_noise(obj, em_read_noise)
            obj.em_read_noise = em_read_noise;
        end

        %% Function to get the EM read-out bias
        % Outputs:
        % read_bias:    Float, EM read-out bias in e-
        function em_read_bias = get.em_read_bias(obj)
            em_read_bias = obj.em_read_bias;
        end

        %% Function to set the EM read-out bias
        % Inputs:
        % read_bias:    Float, EM read-out bias in e-
        function set.em_read_bias(obj, em_read_bias)
            obj.em_read_bias = em_read_bias;
        end

        %% Function to get the rms standard read-out noise
        % Outputs:
        % read_out_noise:   Float, rms standard read-out noise in e-
        function std_read_noise = get.std_read_noise(obj)
            std_read_noise = obj.std_read_noise;
        end

        %% Function to set the rms standard read-out noise
        % Inputs:
        % read_out_noise:   Float, rms standard read-out noise in e-
        function set.std_read_noise(obj, std_read_noise)
            obj.std_read_noise = std_read_noise;
        end

        %% Function to get the standard read-out bias
        % Outputs:
        % read_bias:    Int, standard read-out bias in e-
        function std_read_bias = get.std_read_bias(obj)
            std_read_bias = obj.std_read_bias;
        end

        %% Function to set the standard read-out bias
        % Inputs:
        % read_bias:    Int, read-out bias in e-
        function set.std_read_bias(obj, std_read_bias)
            obj.std_read_bias = std_read_bias;
        end

        %% Function to get the mean gain
        % Outputs:
        % mean_gain:    Float, mean gain of the EM register
        function mean_gain = get.mean_gain(obj)
            mean_gain = obj.mean_gain;
        end

        %% Function to set the mean gain
        % Inputs:
        % mean_gain:    Float, mean gain of the EM register
        function set.mean_gain(obj, mean_gain)
            obj.mean_gain = mean_gain;

            obj.set_em_mult_prob();
        end

        %% Function to define the multiplication probability in the EM
        % register. Called when setting the mean gain of the detector
        function set_em_mult_prob(obj)
            mult_prob = obj.mean_gain.^(1 / obj.num_mult) - 1;

            % Set the multiplication properties for the
            % electron-multiplying elements
            [obj.srl_reg(end - obj.num_mult + 1 : end).mult_prob] = ...
                deal(mult_prob);
        end

        %% Function to retrieve the length of the serial register.
        % Outputs:
        % srl_len:  Int, length of the serial register in pixels
        function srl_len = get_srl_len(obj)
            srl_len = obj.num_std_os + obj.par_wid + obj.num_ad_elem + ...
                obj.num_mult;
        end

        %% Function to retrieve a subsection of the parallel section. Note
        % that the coordinate system origin is located at the bottom-left
        % corner (closest to the standard amplifier)
        % Inputs:
        % x1:   Int, first x-coordinate, must be less than or equal to x2
        %       (1 - img_wid)
        % y1:   Int, first y-coordinate, must be less than or equal to y2
        %       (1 - img_len)
        % x2:   Int, second x-coordinate, must be greater than or equal to 
        %       x1 (1 - img_wid)
        % y2:   Int, second y-coordinate, must be greater than or equal to
        %       y2 (1 - img_wid)
        function pix_array = get_sub_par(obj, x1, y1, x2, y2)
            if x1 <= x2 && y1 <= y2 && x1 >= 1 && y1 >= 1 ...
                    && x2 <= obj.par_wid && y2 <= obj.par_len
                pix_array = obj.par_sec(y1 : y2, x1: x2);
            else
                error(['Detector::get_sub_par: Provided coordinates are'...
                    ' out-of-bounds']);
            end
        end

        %% Function to retrieve a subsection of the serial register. Note
        % that indexing begins at the first pixel before the standard
        % output amplifier
        % Inputs:
        % x1:   Int, first coordinate, must be less than or equal to x2
        %       (1 - num_std_os + img_wid + num_em_os + num_mult)
        % x2:   Int, second x-coordinate, must be greater than or equal to 
        %       x1 (1 - num_std_os + img_wid + num_em_os + num_mult)
        function pix_array = get_sub_srl(obj, x1, x2)
            if x1 <= x2 && x1 >= 1 && x2 <= obj.get_srl_len()
                pix_array = obj.srl_reg(x1 : x2);
            else
                error(['Detector::get_sub_srl: Provided coordinates are'...
                    ' out-of-bounds']);
            end
        end

        %% Function to erase all charge in the parallel and serial register
        % sections
        function clear_detector(obj)
            [obj.par_sec.charge_cloud] = deal(0);
            [obj.srl_reg.charge_cloud] = deal(0);
        end

        %% Function to simulate an integration period on the detector.
        % simulates generation of photoelectrons and dark current
        % Inputs:
        % photons:  Float array, mean photons sec^-1 incident on the pixels
        %           of the detector. Must be the same size as the parallel
        %           section. The third dimension tracks photon maps for
        %           different wavelengths
        % wl:       Float vector, wavelength of the incoming photons in nm
        %           corresponding to the third dimension of photons input.
        % time:     Float, time in seconds over which to generate charge
        %           and dark current signal on the parallel section.
        function integrate(obj, photons, wl, time)
            % Generate charge in the detector for each photon map.
            for i = 1 : length(wl)
                % Only perform charge generation if light is actually
                % falling on the detector
                if any(any(photons(:, :, i)))
                    obj.generate(photons(:, :, i), wl(i), time);
                end
            end

            % Generate dark current over the integration time
            obj.parallel_dark(time);
            obj.serial_dark(time);
        end

        %% Function to generate photoelectrons in the parallel section
        % Inputs:
        % photons:  Float array, mean photons sec^-1 incident on the pixels
        %           of the detector. Must be the same size as the parallel
        %           section.
        % wl:       Float, wavelength of the incoming photons in nm
        % time:     Float, time in seconds over which to generate charge
        %           and dark current signal on the parallel section.
        function generate(obj, photons, wl, time)
            if all(size(photons) == [obj.par_len obj.par_wid])
                arrayfun(@(x, photons) x.generate(photons, wl), ...
                    [obj.par_sec], poissrnd(photons * time));
            else
                error(['Detector::generate: Input photon flux map is '...
                    'incorrect size']);
            end
        end

        %% Function to generate dark current in the parallel section
        % Inputs:
        % time: Float, time in seconds over which to generate dark current.
        function parallel_dark(obj, time)
            arrayfun(@(x) x.dark_gen(time), [obj.par_sec]);
        end

        %% Function to generate dark current in the serial register section
        % Inputs:
        % time: Float, time in seconds over which to generate charge and
        %       dark current signal in the serial register.
        function serial_dark(obj, time)
            arrayfun(@(x) x.dark_gen(time), [obj.srl_reg]);
        end

        %% Function to forward clock the parallel section. Accounts for
        % photon production during this time in both the parallel and
        % serial register
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section.
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % vert_freq:    Float, vertical transfer frequency in Hz
        % opt:          Boolean, controls whether to move charge into the
        %               serial register. Set to true to move charge into
        %               the serial register, otherwise do not move charge
        %               in the first row
        function fw_clock_par(obj, photons, wl, vert_freq, opt)
            obj.integrate(photons, wl, 1 / vert_freq);
            
            if opt
                % Collect charge from each pixel being transferred
                charge_arr = arrayfun(@(x) x.transfer_out(), ...
                    [obj.par_sec]);
    
                % Move charge forward by one pixel
                arrayfun(@(x, c) x.transfer_in(c), ...
                    [obj.par_sec(1 : end - 1, :)], charge_arr(2 : end, :));
    
                % Move charge into serial register
                arrayfun(@(x, c) x.transfer_in(c), ...
                    [obj.srl_reg(obj.num_std_os + 1 : obj.num_std_os + ...
                    obj.par_wid)], charge_arr(1, :));
            else
                % Collect charge from each pixel being transferred
                charge_arr = arrayfun(@(x) x.transfer_out(), ...
                    [obj.par_sec(2 : end, :)]);
    
                % Move charge forward by one pixel
                arrayfun(@(x, c) x.transfer_in(c), ...
                    [obj.par_sec(1 : end - 1, :)], charge_arr);
            end
        end

        %% Function to reverse clock the parallel section. Accounts for
        % photon production during this time in both the parallel and
        % serial register
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section.
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % vert_freq:    Float, vertical transfer frequency in Hz
        % opt:          Boolean, controls whether to move charge out of the
        %               serial register. Set to 1 move out of the serial
        %               register, otherwise charge is not removed from the
        %               serial register
        function rv_clock_par(obj, photons, wl, vert_freq, opt)
            obj.integrate(photons, wl, 1 / vert_freq);
            
            if opt
                % Collect charge from each pixel being transferred
                charge_arr = arrayfun(@(x) x.transfer_out(), ...
                    [obj.par_sec(1 : end - 1, :)]);
                charge_arr_srl = arrayfun(@(x) x.transfer_out(), ...
                    [obj.srl_reg(obj.num_std_os + 1 : ...
                    obj.num_std_os + obj.par_wid)]);

                % Move charge backwards by one pixel
                arrayfun(@(x, c) x.transfer_in(c), ...
                [obj.par_sec(2 : end, :)], charge_arr);
            
                % Move charge out of serial register
                arrayfun(@(x, c) x.transfer_in(c), ...
                    [obj.par_sec(1, :)], charge_arr_srl);
            else
                % Collect charge from each pixel being transferred
                charge_arr = arrayfun(@(x) x.transfer_out(), ...
                    [obj.par_sec(1 : end - 1, :)]);

                % Move charge backwards by one pixel
                arrayfun(@(x, c) x.transfer_in(c), ...
                [obj.par_sec(2 : end, :)], charge_arr);
            end
            
        end

        %% Function to forward clock the serial register (towards the EM
        % output). Note that this function is different from em_read in 
        % that charge is not removed from the last register for read out,
        % but instead stays put. Typically em_read will be used for most
        % applications.
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section. Can be left unspecified so
        %               that the function only generates serial dark
        %               current, and parallel charge generation can be
        %               performed later
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % horz_freq:    Float, horizontal trasnfer frequency in Hz
        function fw_clock_srl(obj, photons, wl, horz_freq)
            if nargin == 3
                obj.integrate(photons, wl, 1 / horz_freq);
            else
                % Must integrate charge generation in the parallel section
                % manually after completing consecutive serial transfers
                obj.serial_dark(1 / horz_freq);
            end

            % Collect charge from each pixel being transferred
            charge_arr = arrayfun(@(x) x.transfer_out(), ...
                [obj.srl_reg(1 : end - 1)]);
            
            % Move charge forward by one pixel
            arrayfun(@(x, c) x.transfer_in(c), [obj.srl_reg(2 : end)], ...
                charge_arr);
        end

        %% Function to reverse clock the serial register (towards the
        % standard output) Note that this function is different from
        % std_read in that charge is not removed from the first register
        % for read out, but instead stays put. Typically std_read will be
        % used for most applications.
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section. Can be left unspecified so
        %               that the function only generates serial dark
        %               current, and parallel charge generation can be
        %               performed later
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % horz_freq:    Float, horizontal trasnfer frequency in Hz
        function rv_clock_srl(obj, photons, wl, horz_freq)
            if nargin == 3
                obj.integrate(photons, wl, 1 / horz_freq);
            else
                % Must integrate charge generation in the parallel section
                % manually after completing consecutive serial transfers
                obj.serial_dark(1 / horz_freq);
            end

            % Collect charge from each pixel being transferred
            charge_arr = arrayfun(@(x) x.transfer_out(), ...
                [obj.srl_reg(2 : end)]);

            % Move charge backward by one pixel
            arrayfun(@(x, c) x.transfer_in(c), ...
                [obj.srl_reg(1 : end - 1)], charge_arr);
        end

        %% Function to forward clock the serial register (towards the EM
        % output) and provide charge for measurement.
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section. Can be left unspecified so
        %               that the function only generates serial dark
        %               current, and parallel charge generation can be
        %               performed later
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % horz_freq:    Float, horizontal trasnfer frequency in Hz
        % Outputs:
        % read_charge:  Double, output from the EM output amplifier in e-
        %               (Note that because the amplifier converts charge to
        %               voltage, this output can have a decimal component).
        function read_charge = em_read(obj, photons, wl, horz_freq)
            if nargin == 3
                obj.integrate(photons, wl, 1 / horz_freq);
            else
                % Must integrate charge generation in the parallel section
                % manually after completing consecutive serial transfers
                obj.serial_dark(1 / horz_freq);
            end

            % Collect charge from each pixel being transferred
            charge_arr = arrayfun(@(x) x.transfer_out(), [obj.srl_reg]);

            % Move charge forward by one pixel
            arrayfun(@(x, c) x.transfer_in(c), [obj.srl_reg(2 : end)], ...
                charge_arr(1 : end - 1));

            % Read the charge stored in the final register, applying 
            % read-out noise and the read-out bias.
            read_charge = double(charge_arr(end)) + ...
                random('Normal', obj.em_read_bias, obj.em_read_noise);
        end

        %% Function to reverse clock the serial register (towards the
        % standard output) and provide charge for measurement from the
        % standard output amplifier.
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section. Can be left unspecified so
        %               that the function only generates serial dark
        %               current, and parallel charge generation can be
        %               performed later
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % horz_freq:    Float, horizontal trasnfer frequency in Hz
        % Outputs:
        % read_charge:  Double, output from the standard output amplifier in
        %               e- (Note that because the amplifier converts
        %               charge to voltage, this output can have a decimal
        %               component).
        function read_charge = std_read(obj, photons, wl, horz_freq)
            if nargin == 3
                obj.integrate(photons, wl, 1 / horz_freq);
            else
                % Must integrate charge generation in the parallel section
                % manually after completing consecutive serial transfers
                obj.serial_dark(1 / horz_freq);
            end

            % Collect charge from each pixel being transferred
            charge_arr = arrayfun(@(x) x.transfer_out(), [obj.srl_reg]);

            % Move charge forward by one pixel
            arrayfun(@(x, c) x.transfer_in(c), ...
                [obj.srl_reg(1 : end - 1)], charge_arr(2 : end));

            % Read the charge stored in the first register, applying
            % read-out noise and read-out bias.
            read_charge = double(charge_arr(1)) + ...
                random('Normal', obj.std_read_bias, obj.std_read_noise);
        end

        %% Function to dump a subsection of the serial register into the
        % dump drain. Note that indexing begins at the first pixel before
        % the standard output amplifier.
        % Inputs:
        % photons:      Float array, mean photons sec^-1 incident on the
        %               pixels of the detector. Must be the same size as
        %               the parallel section.
        % wl:           Float vector, wavelength of the incoming photons in
        %               nm corresponding to the third dimension of photons
        %               input.
        % x1:           Int, first coordinate, must be less than or equal
        %               to x2 (1 - num_std_os + img_wid + num_em_os +
        %               num_mult)
        % x2:           Int, second x-coordinate, must be greater than or
        %               equal to x1 (1 - num_std_os + img_wid + num_em_os +
        %               num_mult)
        % vert_freq:    Float, vertical transfer frequency in Hz
        % opt:          Boolean, set to 1 for a zero time transfer (no dark
        %               current generation, or 0 to set the time to that
        %               specified by the vertical frequency.
        function dump_drain(obj, photons, wl, x1, x2, vert_freq, opt)
            if x1 >= 1 && x2 <= obj.get_srl_len() && x1 <= x2
                if ~opt
                    obj.integrate(photons, wl, 1 / vert_freq);
                end

                % Collect charge from each pixel being transferred
                arrayfun(@(x) x.transfer_out(), obj.srl_reg(x1 : x2));
            else
                error(['Detector::dump_drain: Provided coordinates are '...
                    'out-of-bounds']);
            end
        end

        %% Function to dump the top row of the parallel section into a dump
        % drain.
        % Inputs:
        % vert_freq:    Float, vertical transfer frequency in Hz
        % opt:          Boolean, set to 1 for a zero time transfer (no dark
        %               current generation, or 0 to set the time to that
        %               specified by the vertical frequency.
        function dump_drain_par(obj, vert_freq, opt)
            if ~opt
                obj.integrate(1 / vert_freq);
            end
            
            % Collect charge from each pixel being transferred
            arrayfun(@(x) x.transfer_out(), obj.par_sec(end, :));
        end
        
        %% Function to return pixel array property maps
        % Inputs:
        % prop:         String, property to retrieve
        % Outputs:
        % par_sec_arr:  Array, array of property values for the parallel
        %               section
        % srl_reg_arr:  Array. array of property values for the serial
        %               register
        function [par_sec_arr, srl_reg_arr] = arr_grab(obj, prop)
            switch prop
                case 'char'
                    par_sec_arr = reshape([obj.par_sec.charge_cloud], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.charge_cloud];
                case 'fwc'
                    par_sec_arr = reshape([obj.par_sec.full_well_capacity], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.full_well_capacity];
                case 'cte'
                    par_sec_arr = reshape([obj.par_sec.charge_transfer_eff], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.charge_transfer_eff];
                case 'dark'
                    par_sec_arr = reshape([obj.par_sec.dark_curr_rate], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.dark_curr_rate];
                case 'cic'
                    par_sec_arr = reshape([obj.par_sec.cic_rate], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.cic_rate];
                case 'mult'
                    par_sec_arr = reshape([obj.par_sec.mult_prob], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.mult_prob];
                case 'bloom'
                    par_sec_arr = reshape([obj.par_sec.bloom_stat()], ...
                        obj.par_len, obj.par_wid);
                    srl_reg_arr = [obj.srl_reg.bloom_stat()];
                otherwise
                    error('Detector::arr_grab: Incorrect property string');
            end
        end

        %% Function to set pixel array properties in the parallel section
        % using input maps
        % Inputs:
        % prop:         String, property to set
        % par_sec_arr:  Array, array of property values for the parallel
        %               section
        function par_arr_set(obj, prop, par_sec_arr)
            if ~all(size(par_sec_arr) == [obj.par_len obj.par_wid])
                error(['Detector::par_sec_arr: Property map is '...
                    'incorrect size']);
            end
            
            % Prep property array for distribution
            p = num2cell(par_sec_arr);

            switch prop
                case 'char'
                    [obj.par_sec.charge_cloud] = deal(p{:});
                case 'fwc'
                    [obj.par_sec.full_well_capacity] = deal(p{:});
                case 'cte'
                    [obj.par_sec.charge_transfer_eff] = deal(p{:});
                case 'dark'
                    [obj.par_sec.dark_curr_rate] = deal(p{:});
                case 'cic'
                    [obj.par_sec.cic_rate] = deal(p{:});
                case 'mult'
                    [obj.par_sec.mult_prob] = deal(p{:});
                otherwise
                    error(['Detector::par_arr_set: Incorrect property'...
                        ' string']);
            end
        end

        %% Function to set pixel array properties in the serial register
        % using input maps
        % Inputs:
        % prop:         String, property to set
        % srl_reg_arr:  Array, array of property values for the parallel
        %               section
        function srl_reg_set(obj, prop, srl_reg_arr)
            if ~all(length(srl_reg_arr) == obj.get_srl_len())
                error(['Detector::srl_reg_arr: Property map is '...
                    'incorrect size']);
            end

            % Prep property array for distribution
            p = num2cell(par_sec_arr);

            switch prop
                case 'char'
                    [obj.srl_reg.charge_cloud] = deal(p{:});
                case 'fwc'
                    [obj.srl_reg.full_well_capacity] = deal(p{:});
                case 'cte'
                    [obj.srl_reg.charge_transfer_eff] = deal(p{:});
                case 'dark'
                    [obj.srl_reg.dark_curr_rate] = deal(p{:});
                case 'cic'
                    [obj.srl_reg.cic_rate] = deal(p{:});
                case 'mult'
                    disp(['Detector::srl_reg_set: WARNING: Setting the '...
                        'multiplication probability directly is not ' ...
                        'recommended. The multiplication probability ' ...
                        'may no longer represent the mean gain']);
                    [obj.srl_reg.mult_prob] = deal(p{:});
                otherwise
                    error(['Detector::srl_reg_set: Incorrect property'...
                        ' string']);
            end
        end

        %% Function to check if the detector is set up
        function bool = check_init(obj)
            bool = ~(isempty(obj.par_sec) || isempty(obj.srl_reg));
        end
    end
end