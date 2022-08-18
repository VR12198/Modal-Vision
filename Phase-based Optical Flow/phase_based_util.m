function [u, v] = phase_based_util(img1, img2, horizontal, vertical, gaborBank, sigma)
    % Gaussion smoothing for noise reduction
    img2 = imgaussfilt(img2, sigma);
    % Full field displacement vectors. u: Horizontal, v: Vertical
    u = zeros(size(img2)); v = u;
    
    if horizontal == 1
        [local_mag_h, local_phase_h] = imgaborfilt(img2,gaborBank(:,1));
        % Sobel gradient method for spatial gradient
        image_grad_x = imgradient(local_phase_h,'sobel');
        % Time difference between frames
        image_gradt_x = local_phase_h - img1.ref_local_phase_h;
        % Create a mask to handle Exploding gradient (too close to zero)
        mask = abs(image_grad_x) .* img1.thres_mask_h > 1e-12;
        % Get horizontal displacement
        u(mask) = image_gradt_x(mask)./...
            image_grad_x(mask);
        u = u.*img1.thres_mask_h;
        
        
    end
    if vertical == 1
        [local_mag_v, local_phase_v] = imgaborfilt(img2,gaborBank(:,2));
        % Sobel gradient method for spatial gradient
        image_grad_y = imgradient(local_phase_v,'sobel');
        % Time difference between frames
        image_gradt_y = local_phase_v - img1.ref_local_phase_v;
        % Create a mask to handle Exploding gradient (too close to zero)
        mask = abs(image_grad_y).* img1.thres_mask_v > 1e-3;
        % Get vertical displacement
        v(mask) = image_gradt_y(mask)./...
            image_grad_y(mask);
        v = v.*img1.thres_mask_v;
    end   
end