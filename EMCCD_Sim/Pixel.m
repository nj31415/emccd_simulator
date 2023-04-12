%% Nicholas Jones - njones31@vt.edu
% The Pixel class is used to represent a single pixel of a detector. Each
% pixel tracks its own charge and supports indepently setting various
% properties related to pixel performance.

classdef Pixel < handle
    properties (Constant)
        e_pair = 3.65;      % Average energy required to produce an e-h
                            % pair in Si at room temperature (eV)
        fano_fac = 0.1;     % Fano factor for Silicon.
    end
    properties (SetAccess = public)
        %% Define pixel properties
        full_well_capacity
        charge_transfer_eff
        dark_curr_rate
        cic_rate
        mult_prob

        charge_cloud = 0;   % Number of electrons stored in the pixel

        % Interacting Quantum Efficiency Tables
        % Can specify interacting quantum efficiency by wavelength. When
        % using charge generation functions, a specified wavelength will be
        % used to interpolate an interacting quantum efficiency for that
        % wavelength. Specified in the range between 0 and 1. Note that
        % this is not the same as quantum efficiency, which is the product
        % of iQE and the quantum yield.

        wl_ref_iqe          % Stores wavelength reference table in nm for
                            % iQE look-up and interpolation
        i_quantum_eff       % Stores iQE corresponding to wavelengths
                            % specified in wl_ref_iqe. Defined by the user
                            % to reflect the device and coatings used

        iqe_intrp_m = 'linear'; % Interpolation method for determining iQE.
                                % Recommend 'linear' or 'nearest'.
    end

    properties (SetAccess = private)
        % Quantum Yield Tables. Implemented based on discussions with
        % NASA JPL collaborators, based on data from 'Characterization of
        % photodiodes as transfer detector standards in the 120 nm to 600
        % nm spectral range' by Kuschnerus et al (1998). At wavelengths
        % greater than 350 nm, quantum yield is assumed to be 1. At
        % wavelengths longer than 123.9 nm, quantum yield is estimated
        % based on equations in Janesick 2001 and McLean 2008. The quantum
        % yield value for 123.9 nm was appended to the data provided by
        % NASA JPL using the equations in Janesick 2001 for quantum yield
        % for photons > 10 eV.

        % Stores wavelength reference table in nm for quantum yield look-up
        % and interpolation
        wl_ref_qy_le = [350 340 330 320 310 300 290 280 270 260 253.7 ...
            238.5 229.6 221.4 218 216 214 213.8 212 210 208 206 204 202 ...
            200 193.7 182.3 175 170 123.9];
        
        % Stores quantum yield corresponding to wavelengths specified in
        % wl_ref_qy_le. Required for photons with energy < 10 eV 
        % (wavelength > 123.9 nm).
        quantum_yield_le = [1 1.004519774 1.009786408 1.019548387 ...
            1.045251748 1.076514563 1.102121951 1.09178125 1.095309735 ...
            1.126615385 1.150846154 1.214636364 1.253649351 1.286662338 ...
            1.300350649 1.308589552 1.320380597 1.321559701 1.332171642 ...
            1.343962687 1.355753731 1.367544776 1.379335821 1.397950495 ...
            1.429237624 1.527792079 1.705214286 1.786035714 1.845480392 ...
            2.739726027];
    end
    
    methods
        %% Consturctor for a Pixel object. When called with no parameters,
        % creates uninstantiated pixels (Needed for array initialization in
        % the detector class). The wl_ref_iqe and i_quantum_eff variables
        % must be inside of a 1x1 cell array in order to make the deal
        % function used in Detector class work properly.
        % Inputs:
        % full_well_capacity:   e- - Full well capacity
        % charge_transfer_eff:  Charge transfer efficiency (0 - 1)
        % cic_rate:             e- pix^-1 tr^-1 - CIC generation rate
        % mult_prob:            Probability of multiplication (0 - 1)
        %                       (Intended for gain registers, for
        %                       non-multiplying pixels set to 0)
        % wl_ref_iqe:           Float vector, reference wavelengths in nm
        %                       for iQE. Must be passed as a vector inside
        %                       a cell
        % i_quantum_eff:        Float vector, interacting quantum
        %                       efficiency. Corresponds to wavelengths
        %                       stored in wl_ref_iqe. Must be passed as a
        %                       vector inside a cell.
        function obj = Pixel(full_well_capacity, charge_transfer_eff, ...
                dark_curr_rate, cic_rate, mult_prob, wl_ref_iqe, ...
                i_quantum_eff)
            if nargin > 0
                if length(wl_ref_iqe) == length(i_quantum_eff)
                    obj.full_well_capacity = full_well_capacity;
                    obj.charge_transfer_eff = charge_transfer_eff;
                    obj.dark_curr_rate = dark_curr_rate;
                    obj.cic_rate = cic_rate;
                    obj.mult_prob = mult_prob;
                    obj.wl_ref_iqe = wl_ref_iqe;
                    obj.i_quantum_eff = i_quantum_eff;
                else
                    error(['Pixel::Pixel:iQE and reference wavelength ' ...
                        'vectors are not the same length.']);
                end
            end
        end

        %% Get the interacting quantum efficiency
        % Outputs:
        % i_quantum_eff : Float vector, interacting quantum efficiency.
        %                 Corresponds to wavelengths stored in wl_ref_iqe.
        function i_quantum_eff = get.i_quantum_eff(obj)
            i_quantum_eff = obj.i_quantum_eff;
        end

        %% Set the interacting quantum efficiency
        % Inputs:
        % i_quauntum_eff:   Interacting quantum efficiency (0 - 1).
        %                   Note that this is not the same as quantum
        %                   efficiency. Corresponds to wavelengths stored
        %                   in wl_ref_iqe.
        function set.i_quantum_eff(obj, i_quantum_eff)
            i_quantum_eff = i_quantum_eff{1};
            if all(i_quantum_eff >= 0 & i_quantum_eff <= 1)
                if isempty(obj.wl_ref_iqe) || ...
                        length(obj.wl_ref_iqe) == length(i_quantum_eff)
                    obj.i_quantum_eff = i_quantum_eff;
                else
                    error(['Pixel::set_i_quantum_eff: Incorrect size of'...
                        ' vector']);
                end
            else
                error(['Pixel::set.i_quantum_eff: Interacting quantum '...
                    'efficiency is out of range (0 - 1)']);
            end
        end

        %% Get the reference wavelengths for the interacting quantum
        % efficiency
        % Outputs:
        % wl_ref_iqe : Float vector, reference wavelengths in nm
        function wl_ref_iqe = get.wl_ref_iqe(obj)
            wl_ref_iqe = obj.wl_ref_iqe;
        end

        %% Set the reference wavelengths for the interacting quantum
        % efficiency
        % Inputs:
        % wl_ref_iqe:   Float vector, reference wavelengths in nm.
        function set.wl_ref_iqe(obj, wl_ref_iqe)
            wl_ref_iqe = wl_ref_iqe{1};
            if isempty(obj.i_quantum_eff) || ...
                    length(wl_ref_iqe) == length(obj.i_quantum_eff)
                obj.wl_ref_iqe = wl_ref_iqe;
            else
                error('Pixel::set.wl_ref_iqe: Incorrect size of vector');
            end
        end

        %% Function to get the interpolation method for the interacting
        % quantum efficiency
        % Outpus:
        % iqe_intrp_m:  String, interpolation method
        function iqe_intrp_m = get.iqe_intrp_m(obj)
            iqe_intrp_m = obj.iqe_intrp_m;
        end

        %% Function to set the interpolation method for the interacting
        % quantum efficiency
        % Inputs:
        % iqe_intrp_m:  String, interpolation method. Must correspond to an
        %               interpolation method for MATLAB's interp1 function.
        function set.iqe_intrp_m(obj, iqe_intrp_m)
            obj.iqe_intrp_m = iqe_intrp_m;
        end

        %% Get the quantum yield data for low energy photons (< 10 eV)
        % Outputs:
        % quantum_yield_le:    Float, quantum yield
        function quantum_yield_le = get.quantum_yield_le(obj)
            quantum_yield_le = obj.quantum_yield_le;
        end

        %% Set the quantum yield for low energy photons (< 10 eV).
        % Inputs:
        % quantum_yield:    Quantum yield ( >= 1, dependent on energy
        %                   of incoming photons)
        function set.quantum_yield_le(obj, quantum_yield_le)
            if all(quantum_yield_le >= 1)
                if isempty(obj.wl_ref_qy_le) || ...
                        length(obj.wl_ref_qy_le) == ...
                        length(quantum_yield_le)
                    obj.quantum_yield_le = quantum_yield_le;
                else
                    error(['Pixel::set.quantum_yield_le: Incorrect ' ...
                        'vector size']);
                end
            else
                error(['Pixel::set.quantum_yield_le: Quantum yield is ' ...
                    'out of range ( < 1)']);
            end
        end

        %% Get the reference wavelengths for the quantum yield data for low
        % energy photons (< 10 eV)
        % Outputs:
        % wl_ref_qy_le:    Float, reference wavelengths in nm
        function wl_ref_qy_le = get.wl_ref_qy_le(obj)
            wl_ref_qy_le = obj.wl_ref_qy_le;
        end

        %% Set the reference wavelengths for the quantum yield for low
        % energy photons (< 10 eV).
        % Inputs:
        % wl_ref_qy_le: reference wavelengths in nm
        function set.wl_ref_qy_le(obj, wl_ref_qy_le)
            if all(wl_ref_qy_le >= 1)
                if isempty(obj.quantum_yield_le) || ...
                        length(obj.quantum_yield_le) == ...
                        length(wl_ref_qy_le)
                    obj.wl_ref_qy_le = wl_ref_qy_le;
                else
                    error(['Pixel::set.wl_ref_qy_le: Incorrect vector ' ...
                        'size']);
                end
            else
                error(['Pixel::set.quantum_yield: Quantum yield is out '...
                    'of range ( < 1)']);
            end
        end
        
        %% Get the full well capacity
        % Outputs:
        % full_well_capacity:   Float, full well capacity in e-
        function full_well_capacity = get.full_well_capacity(obj)
            full_well_capacity = obj.full_well_capacity;
        end

        %% Set the full well capacity
        % Inputs:
        % full_well_capacity:   uint32, full well capacity in e-
        function set.full_well_capacity(obj, full_well_capacity)
            if full_well_capacity >= 0
                obj.full_well_capacity = uint32(full_well_capacity);
            else
                error(['Pixel::set.full_well_capacity: Full well '...
                    'capacity is out of range (> 0)']);
            end
        end

        %% Get the charge transfer efficiency
        % Outputs:
        % charge_transfer_eff:  Float, charge transfer efficiency
        function charge_transfer_eff = get.charge_transfer_eff(obj)
            charge_transfer_eff = obj.charge_transfer_eff;
        end

        %% Set the charge transfer efficiency
        % Inputs:
        % charge_transfer_eff:  Charge transfer efficiency (0 - 1)
        function set.charge_transfer_eff(obj, charge_transfer_eff)
            if charge_transfer_eff >= 0 && charge_transfer_eff <= 1
                obj.charge_transfer_eff = charge_transfer_eff;
            else
                error(['Pixel::set.charge_transfer_eff: Charge transfer'...
                    ' is out of range (0 - 1)']);
            end
        end

        %% Get the dark current generation rate
        % Outputs:
        % dark_curr_rate:  Float, dark current rate in e- pix^-1 sec^-1
        function dark_curr_rate = get.dark_curr_rate(obj)
            dark_curr_rate = obj.dark_curr_rate;
        end

        %% Set the dark current generation rate
        % Inputs:
        % dark_curr_rate:  Float, dark current rate in e- pix^-1 sec^-1
        function set.dark_curr_rate(obj, dark_curr_rate)
            obj.dark_curr_rate = dark_curr_rate;
        end

        %% Get the clock-induced-charge generation rate
        % Outputs:
        % cic_rate:  Float, clock-induced-charge rate in e- pix^-1 tr^-1
        function cic_rate = get.cic_rate(obj)
            cic_rate = obj.cic_rate;
        end

        %% Set the clock-induced-charge generation rate. Note this is in
        % terms of per transfer (between pixels), rather than per frame
        % Inputs:
        % cic_rate:  Float, dark current rate in e- pix^-1 tr^-1
        function set.cic_rate(obj, cic_rate)
            obj.cic_rate = cic_rate;
        end

        %% Get the multiplication probability
        % Outputs:
        % mult_prob:  Float, multiplication probability
        function mult_prob = get.mult_prob(obj)
            mult_prob = obj.mult_prob;
        end

        %% Set the multiplication probability. This is used to set the
        % probability of multiplication in serial registers
        % Inputs:
        % mult_prob:            Probability of multiplication (0 - 1)
        %                       (Intended for gain registers, for
        %                       non-multiplying pixels set to 0)
        function set.mult_prob(obj, mult_prob)
            if mult_prob >= 0 && mult_prob <= 1
                obj.mult_prob = mult_prob;
            else
               error(['Pixel::set.mult_prob: multiplication probability'...
                   ' is out of range (0 - 1)']);
            end
        end

        %% Get the current size of the charge cloud
        % Outputs:
        % charge_cloud: Int, size of chare cloud in e-
        function charge_cloud = get.charge_cloud(obj)
            charge_cloud = obj.charge_cloud;
        end

        % Set the size of the charge cloud
        % Inputs:
        % charge_cloud: Int, size of chare cloud in e-
        function set.charge_cloud(obj, charge_cloud)
            obj.charge_cloud = charge_cloud;

            obj.clear_excess();
        end

        %% Function to get blooming state of the pixel
        % Outputs:
        % bloom:    boolean, true if the amount of charge stored in the
        % pixel exceeds the full well capacity, false otherwise
        function bloom_stat = bloom_stat(obj)
            bloom_stat = obj.charge_cloud > obj.full_well_capacity;
        end
        
        %% Function to generate charge based on incoming photons.
        % Calculates the number of electrons generated by taking into
        % account the interacting quantum efficiency, quantum yield, and
        % Fano factor. Should only be called on image section pixels.
        % Inputs:
        % photons:  Integer, number of photons incident on the pixel
        % wl:       Float, wavelength of the incoming photons in nm
        function generate(obj, photons, wl)
            % Determine the quantum efficiency to use. If the wavelength
            % passed in is outside the range specified in wl_ref_iqe,
            % defaults to 0.
            qe_intrp = interp1(obj.wl_ref_iqe, obj.i_quantum_eff, wl, ...
                obj.iqe_intrp_m, 0);
            % Determine which photons interact
            interact = obj.binom_rnd(photons, qe_intrp);
            
            % Only continue if there are interacting photons. Otherwise do
            % nothing.
            if interact > 0
                % Determine the quantum yield to use and calculate how much
                % charge to add to the charge cloud. The form of the 
                % calculation depends on the energy of the photon. 
                % Also, adds electrons produced by
                % interacting photons to the charge cloud, as well as apply
                % Fano noise according to model described in 'NUV
                % Performance of e2v Large BICMOS Array for CASTOR' (Scott 
                % et al, 2016) for low energy (< 10 eV) photons. Allows for
                % generation of up to 2 additional electrons.
                if wl <= 123.9
                    % Calculate the mean quantum yield based on equation
                    % 1.2 from "Scientific Charge-Coupled Devices" by
                    % Janesick, 2001 for photons with sufficiently high
                    % energy (> 10 eV, or equivalently wavelength > 123.9
                    % nm)
                    ph_nrg = obj.wl_2_nrg(wl);
                    qy_intrp = ph_nrg / obj.e_pair;

                    % Determine the standard deviation of the spread in
                    % electron - hole pairs generated, in e-. Uses equation
                    % 12.6 from "Electronic Imaging in Astronomy" by
                    % McLean, 2008. Note that this calculates the standard
                    % deviation rather than the FWHM as given in the
                    % source.
                    sigma_e = sqrt(obj.fano_fac * ph_nrg / obj.e_pair);

                    % All interacting photons generate 1 electron - hole
                    % pair, so the mean used in the gaussian distribution
                    % is 1 less than the quantum yield. Round function is
                    % used to remove fractional component of generated
                    % numbers while still approximating a Gaussian
                    % distribution.
                    qy_e = round(normrnd(qy_intrp - 1, sigma_e, ...
                        interact, 1));
                    
                    % Removes potential negatives generated from the normal
                    % distribution, though this should be extremely
                    % unlikely for > 10 eV photons
                    qy_e(qy_e < 0) = 0;
                    
                    % Add electrons to the charge cloud
                    obj.charge_cloud = obj.charge_cloud + interact + ...
                        sum(qy_e);
                elseif wl > max(obj.wl_ref_qy_le)
                    % Quantum yield is 1 for sufficiently low energy
                    % photons. Just add the number of interacting photons
                    % to the charge cloud.
                    obj.charge_cloud = obj.charge_cloud + interact;
                else
                    % Interpolate the quantum yield from the reference data
                    % for intermediate energy photons (where quantum yield
                    % is greater than 1, but energy is less than 10 eV).
                    qy_intrp = ...
                        interp1(obj.wl_ref_qy_le, obj.quantum_yield_le, ...
                        wl, 'linear');

                    % Calculate probability of generating an additional e-h
                    % pair, following the procedure laid out in 'NUV
                    % Performance of e2v Large BICMOS Array for CASTOR'
                    % (Scott et al, 2016).
                    p_prob = -0.5 + 0.5 * sqrt(-3 + 4 * qy_intrp);
                    
                    % Allows the photon to attempt to generate up to 2
                    % additional e-h pairs
                    obj.charge_cloud = obj.charge_cloud + interact + ...
                        obj.binom_rnd(interact, p_prob) + ...
                        obj.binom_rnd(interact, p_prob^2);
                end

                obj.clear_excess();
            end
        end

        %% Function to generate dark current and add it to the charge cloud
        % Inputs:
        % time: time over which dark current is generated in seconds
        function dark_gen(obj, time)
            % Like CIC, dark current is modeled as a Poisson process with
            % average value given by the time multiplied by the dark
            % current generation rate.
            obj.charge_cloud = obj.charge_cloud + ...
                poissrnd(time * obj.dark_curr_rate);

            obj.clear_excess();
        end

        %% Function to transfer charge from the charge cloud. Calculates
        % the number of electrons to transfer based on the CTE, removes
        % those electrons from the charge cloud before returnnig them as
        % part of this function.
        % Output:
        % charge:   Integer, number of electrons transfered out of the
        %           pixel and removed from the charge cloud
        function charge = transfer_out(obj)
            charge = obj.binom_rnd(obj.charge_cloud, ...
                obj.charge_transfer_eff);
            obj.charge_cloud = obj.charge_cloud - charge;

            if obj.charge_cloud < 0
                error('Pixel::transfer: Charge cloud has become negative');
            end
            
            charge = charge + obj.cic_gen();
        end
        
        %% Function to generate CIC. CIC is modelled as a Poisson process 
        % where the mean average rate is given by the measured detector CIC
        % rate. Note: From literature, it seems that most CIC originates in
        % the serial register, specifically the EM register, so it may be 
        % approriate to set a variable CIC over the detector area to
        % reflect this.
        % Outputs:
        % cic:  Int, number of CIC e- generated
        function cic = cic_gen(obj)
            cic = poissrnd(obj.cic_rate);
        end

        %% Function to transfer charge into the charge cloud. If in the
        % serial register, electron multiplication also occurs here.
        % Inputs:
        % charge:   Integer, number of electrons transferred into the pixel
        %           and added to the charge cloud.
        function transfer_in(obj, charge)
            obj.charge_cloud = obj.charge_cloud + charge;
            
            % Perform charge multiplication
            if obj.mult_prob > 0
                obj.em_mult();
            end

            obj.clear_excess();
        end

        %% Function to calculate electron multiplication on the current
        % charge cloud present in the pixel. Adds the multiplied charge to
        % the charge cloud
        function em_mult(obj)
            obj.charge_cloud = obj.charge_cloud + ...
                obj.binom_rnd(obj.charge_cloud, obj.mult_prob);
            obj.clear_excess();
        end

        %% Function to inject charge into the pixel
        % Inputs:
        % inj_charge:   Charge to inject in the pixel, in e-
        function inject(obj, inj_charge)
            obj.charge_cloud = obj.charge_cloud + inj_charge;
            obj.clear_excess();
        end

        %% Function to clear excess charge during pixel blooming conditions.
        % Note that any excess charge is removed from the pixel when this
        % function is called. This function is called at the end of any
        % process that generates charge in the pixel to roughly simulate
        % saturation.
        % Outputs:
        % excess:   Integer, excess charge stored in the pixel beyond full
        %           well capacity
        function excess = clear_excess(obj)
            if obj.charge_cloud > obj.full_well_capacity
                excess = obj.charge_cloud - obj.full_well_capacity;
                obj.charge_cloud = obj.full_well_capacity;
            else
                excess = 0;
            end
        end

        %% Function to generate binomial random numbers. Faster method than
        % binornd from the Statistics Toolbox. Implements gpuArray for
        % number of trials greater than 20000 (GPU based processing becomes
        % faster at this point according to some basic comparison testing
        % using timeit).
        % Inputs:
        % n:    Integer - number of trials
        % p:    Float - probability of success
        % Outputs:
        % rnd:  Integer - random number from the binomial distribution
        function rnd = binom_rnd(~, n, p)
            if n > 20000
                rnd = gather(sum(gpuArray.rand(n, 1) < p));
            else
                rnd = sum(rand(n, 1) < p);
            end
        end

        %% Function to convert between photon wavelength and energy in eV.
        % Conversion based on Eq. 1.1 from "Scientific Charge-Coupled
        % Devices" by Janesick, 2001
        % Inputs:
        % lambda:   Float - photon wavelength in nm
        % Outputs:
        % energy:   Float - photon energy in eV
        function energy = wl_2_nrg(~, lambda)
            energy = 1239 / lambda;
        end

        %% Function to convert between photon energy in eV and wavelength
        % in nm. Conversion based on Eq. 1.1 from "Scientific
        % Charge-Coupled Devices" by Janesick, 2001
        function wl = nrg_2_wl(~, nrg)
            wl = 1239 / nrg;
        end
    end
end