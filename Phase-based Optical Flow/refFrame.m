classdef refFrame
    properties
        image
        thres_mask_h
        thres_mask_v
        ref_local_phase_h
        ref_local_phase_v
    end
    methods
        function obj = refFrame(img, gaborBank, sigma, dilation)
            if isempty(dilation)
                dilation = 0;
            end
            obj.image = imgaussfilt(img, sigma);
            
            % Get magnitude and phase of gabor filter output.
            % Output is a vector of 2D matrix (frame) for each orientation
            % in the filter bank.
            [reference_mag, reference_phase] = imgaborfilt(img, gaborBank);
            
            k = 60; % median k highest value of reference_mag used for threshold calculation
            % Dilation kernel 
            se = strel('square',2);
            
            [x,y] = size(reference_mag(:,:,1));
            % Threshold below which mask will be 0
            thres_h = median(maxk(reshape(reference_mag(10:x-10,10:y-10,1),[],1),k));
            obj.thres_mask_h = ones([x,y]);
            obj.thres_mask_h(reference_mag(:,:,1) < thres_h) = 0;
            if dilation == 1
                obj.thres_mask_h = imdilate(obj.thres_mask_h, se);
            end
            
            % Threshold below which mask will be 0
            thres_v = median(maxk(reshape(reference_mag(10:x-10,10:y-10,2),[],1),k));
            obj.thres_mask_v= ones([x,y]);
            obj.thres_mask_v(reference_mag(:,:,2) < thres_v) = 0;
            if dilation == 1
                obj.thres_mask_v = imdilate(obj.thres_mask_v, se);
            end

            obj.ref_local_phase_h = reference_phase(:,:,1);
            obj.ref_local_phase_v = reference_phase(:,:,2);
        end
    end
end
            
            
