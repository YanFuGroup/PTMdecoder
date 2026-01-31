classdef CReviewPSM
%CREVIEWPSM review the peptide(s)-spectrum match
%   Show how good the match is.
%   Plot the the peaks and annotate the matched peak.
%   Give a weighted score to help the researcher finding good matches.
    
    properties
        m_peptides;     % struct: [seq, mod_mass, mod_pos]
%               seq:    the sequence of the peptide
%               mod_mass:   an array of modification names,
%                   e.g. [42.010565,42.010565,42.010565,56.026215]
%               mod_pos:    an array of modification position on peptide
%                   e.g. [2,5,9,13]
        m_spectrum;     % struct: [peaks, pre_charge];
%               peaks:      [m/z1, intens1; m/z2, intens2;...]
%               pre_charge: the precursor charge according to the spectrum
        m_tolerance;    % struct: [value, is_ppm]
%               value:      the value of mass tolerance
%               is_ppm:     if the type of mass tolerance is ppm, boolean
        m_main_score;   % main score provided by search engine
        m_is_high_score_better;
        m_main_score_upper_bound;
        m_main_score_lower_bound;
        m_weight_factor;% different weight for five different score

        m_all_match_ions;   % double: [match_type, match_pos, charge, expe_which, pep_which]
%               match_type: the length is equal to the number of matched 
%                   peaks, 1 means b-ion, 2 means y-ion
%               match_pos:  same length as match_type, each element shows 
%                   the position of the ion on sequence
%               charge:     the charge of theoretical ion
%               expe_which: the index of matched experimental peak in 
%                   experimental spectrum
%               pep_which:  the index of matched peptide in peptides list
    end
    
    methods
        function obj = CReviewPSM(varargin)
            % Input (3 args):
            %   peptides (1 x P struct)
            %       peptide(s) with fields: seq (char), mod_mass (1 x M double), mod_pos (1 x M double)
            %   spectrum (struct)
            %       fields: peaks (N x 2 double [m/z, intensity]), pre_charge (1 x 1 double/int), pre_mz (1 x 1 double)
            %   tolerance (struct)
            %       fields: value (double), is_ppm (logical)
            % Input (8 args):
            %   peptides, spectrum, tolerance, main_score, is_high_score_better,
            %   main_score_upper_bound, main_score_lower_bound, weight_factor (1 x 5 double)
            if nargin ~= 3 && nargin ~= 8
                error('Wrong number of input arguments. 3 or 8 arguments are needed.');
            end
            if nargin == 3
                obj.m_peptides = varargin{1};
                obj.m_spectrum = varargin{2};
                obj.m_tolerance = varargin{3};
            elseif nargin == 8
                obj.m_peptides = varargin{1};
                obj.m_spectrum = varargin{2};
                obj.m_tolerance = varargin{3};
                obj.m_main_score = varargin{4};
                obj.m_is_high_score_better = varargin{5};
                obj.m_main_score_upper_bound = varargin{6};
                obj.m_main_score_lower_bound = varargin{7};
                weight_factor = varargin{8};
                if ~isequal([1,5], size(weight_factor))
                    error('The weight factors are not in [1*5] format.');
                end
                if sum(weight_factor) == 0
                    error('The sum of weight factors cannot be 0.');
                end
                weight_factor = weight_factor/sum(weight_factor);
                obj.m_weight_factor = weight_factor;
            end

            filter_alpha = 0.01;
            filter_thres = filter_alpha*max(obj.m_spectrum.peaks(:,2));
            obj.m_spectrum.peaks = obj.m_spectrum.peaks( ...
                obj.m_spectrum.peaks(:,2)>filter_thres,:);
            obj = match_allp_1s(obj);
        end
        
        % Match all theoretical peptide ions with peaks in the spectrum.
        obj = match_allp_1s(obj);

        % Match one theoretical peptide ions with peaks in the spectrum.
        [match_ions] = match_1p_1s(obj, pep_seq, spec_peaks, tolerance);

        % Give a weighted score to help the researcher finding good matches.
        review_score = reviewScoreing(obj);

        % Plot the spectra and mark the matched peaks in the spectrum
        plotMatch(obj);

        
    end
end

