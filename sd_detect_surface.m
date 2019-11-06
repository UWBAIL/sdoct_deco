function surface_z = sd_detect_surface(img, img_xz, windowlength, maxjump, minseg)

    % Construct Derivative-of-Gaussian operator
    dog = diff(gausswin(windowlength + 1));

    % Storage for outputs
    [nz, nt, nx] = size(img);


    % Preprocessing - median filter in axial direction
    for j = 1:nx 
        filtered_line = medfilt1(img_xz(:,j), 5);
        img_xz(3:(nz-2),j) = filtered_line(3:(nz-2));
    end
    
    % Apply DoG differentiation to each column of the image array,
    % Locate the surface as the maximum derivative point. Because 
    % the convolution is not valid near the edges, shift the surface 
    % indices to account for the DoG window length.
    img_deriv = conv2(dog, 1, img_xz, 'valid');
    [d, sidx] = max(img_deriv, [], 1);
    surface_z = sidx + windowlength;

    
%% The following section applies to complicated/curved surfaces (in imaging mode- with a FOV in x.
%% Since this code is at a single location, the surface detection is not necessary

%     % Apply median filter to remove some outlier points from surface
%     % identification. After median filtering the corner points may be 
%     % incorrect as a result of zero padding - we therefore ignore the
%     % filtered corner points.
%     zfilt = medfilt1(surface_z, 5);
%     surface_z(3:nx-2) = zfilt(3:nx-2);
% 
%     % Perform Gaussian filtering to smooth the surface mesh
%     for i = 1:5
%         surface_z = imgaussfilt(surface_z);
%     end
% 
%     % Reject jumps in the surface larger than maxjump
%     zx = abs(diff(surface_z));
%     surface_z(zx > maxjump) = NaN;
% 
%     % Search the surface and discard contiguous segments that are
%     % smaller than minseg.
%     idx = 1;
%     while idx <= nx
%         % Initialize segment
%         segment_start = idx;
%         segment_length = 1;
% 
%         % Scan until we find a NaN, keeping track of the segment length
%         while (idx <= nx) && ~isnan(surface_z(idx))
%             segment_length = segment_length + 1;
%             idx = idx + 1;
%         end
% 
%         % We've found a NaN or reached the end, test segment length
%         if segment_length < minseg
%             surface_z(segment_start:(idx-1)) = NaN;
%         end
% 
%         % Scan until we find another valid non-NaN segment
%         while (idx <= nx) && isnan(surface_z(idx))
%             idx = idx + 1;
%         end
%     end
% 
%     % Fill remaining gaps within the support of the data
%     first_valid = find(~isnan(surface_z), 1, 'first');
%     last_valid = find(~isnan(surface_z), 1, 'last');
%     warning('off', 'MATLAB:interp1:NaNstrip'); 
%     surface_z(first_valid:last_valid) = interp1(surface_z, first_valid:last_valid, 'spline');
%     warning('on', 'MATLAB:interp1:NaNstrip'); 
