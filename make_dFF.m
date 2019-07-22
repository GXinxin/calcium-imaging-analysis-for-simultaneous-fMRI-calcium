function make_dFF(data_dir, results_dir, fileInfo, makeTiffMovie)

% Function to make deltaF/F pixel-by-time array after detrending data by
% tophat filtering and regressing Ca-independent (violet excited)
% signal from the Ca-dependent (blue excited) signal.
%
% Inputs:
%     tiff_base: the base file name for the raw widefield tiff files
%     blueInitial: the index of the first frame that was collected with blue excitation
%     mask_fn: the filename of the mask (ImageJ .roi file) to remove all of the non-brain parts of the tiff files
%     rotateAngle: the angle to rotate the image so that the central sinus is vertical
%



% initialize parameters
isTranslate2Initial = 1; % register each frame to the second frame of all the recordings
waveLengthId = [3]; % 1 - blue only, 2 - UV only, 3 - regressed blue
waveLengthfn = {'blue', 'UV', 'regress'};
hat = 300; % window size for rolling hat algorithm
isTranslate = 1; % whether to motion correct the data by rigid translation

% get the number of movies from fileInf
no_movies = size(fileInfo, 1)



for n = 1 : no_movies
    
    fprintf('Preprocess movie %d. \n', n)
    tic
    clear mask_id filtered1 filtered2 A_filtered B_filtered diff
    
    % load movie information from fileInfo
    fn = fileInfo{n, 1};
    blueInitial = fileInfo{n, 2}; % 1: odd frames are blue channels frames, 2: even frames
    mask_fn = fileInfo{n, 3}; % a hand-drawn binary brain mask
    rotateAngle = fileInfo{n, 4}; % rotate the movie so it is upright
    downSampleRatio = fileInfo{n, 10}; % down sample ratio for preprocessing
    range_r = cell2mat(fileInfo(n, 5:6)); % row range of an ROI including the fluorescent bead for motion correction
    range_c = cell2mat(fileInfo(n, 7:8)); % column range of an ROI including the fluorescent bead for motion correction
    range_r = ceil(range_r * downSampleRatio); % bead ROI multiplied by the corresponding down sample ratio
    range_c = ceil(range_c * downSampleRatio);
    isStim = fileInfo{n, 9}; % whether has stimulation
    if isStim
        stim_fnm = fileInfo{n, 11}; % the spike2 filename with stimulation profiles
    end
    x_res = 0.025 / downSampleRatio; % original pixel resolution = 25*25um^2
    y_res = 0.025 / downSampleRatio;
        
    
    
    % Open .tiff files and concatenate separate arrays for blue and violet
    % frames.
    concatList = dir(fullfile([data_dir, '/', fn]));
    
    A = [];
    B = [];
    for c = 1 : length(concatList)
        fprintf('Opening movie segment %d out of %d. \n',c,length(concatList))
        imgall = openMovie([concatList(c).folder,'/',concatList(c).name]);
        szall = size(imgall);
        if mod(c, 2) == 1
            blueFrames = blueInitial : 2: szall(3);
            uvFrames = setdiff(1:szall(3), blueFrames);
        else
            uvFrames = blueInitial : 2: szall(3);
            blueFrames = setdiff(1:szall(3), uvFrames);
        end
        
        A = cat(3, A, imresize(imgall(:, :, blueFrames), downSampleRatio, 'bilinear'));
        B = cat(3, B, imresize(imgall(:, :, uvFrames), downSampleRatio, 'bilinear'));
        
        
        wrongId = 200:200:size(B, 3); % correct frames at each scanner trigger
        for id = wrongId(1:end-1)
            A(:, :, id) = (A(:, :, id-1) + A(:, :, id+1))/2;
            B(:, :, id) = (B(:, :, id-1) + B(:, :, id+1))/2;
        end
        
        % use a subset of frames to generate an average frame for registration
        if n * c == 1
            data4avg = imgall(:, :, blueFrames);
        end
        clear imgall
    end
    
    
    
    % generate a raw data movie for demonstration
    if makeTiffMovie
        fprintf('generating a raw movie for blue channel...')
        tmp = imrotate(A(:, :, 1:2087), 90); % the movie in the supplementary figure was rotated for visualization purpose
        tmp = 255 * mat2gray(tmp) + 1;
        frStart = 1;
        frEnd = size(tmp, 3);
        Iarr2avi_grayScale(tmp, frStart, frEnd, [results_dir, '/', fn(1:end-5), '_submovie.avi']); % dF/F with stim marker
        clear tmp
    end
    
    
    
    % match number of frames in A and B
    if size(A, 3) >= size(B, 3)
        A = A(:, :, 1:size(B, 3));
    else
        B = B(:, :, 1:size(A, 3));
    end
    
    
    
    % translate individual frames based on
    fprintf('running motion correction')
    if isTranslate
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumIterations = 100;
        
        
        % blue channel
        if n == 1
            fixed_A = A(:, :, 2);
        else
            if ~isTranslate2Initial
                fixed_A = A(:, :, 2);
            end
        end
        subMovie = A(range_r(1):range_r(2), range_c(1):range_c(2), :);
        sub_fixed = fixed_A(range_r(1):range_r(2), range_c(1):range_c(2));
        for f = 2 : size(A, 3)
            tform = imregtform(sub_fixed, subMovie(:, :, f), ...
                'translation', optimizer, metric);
            tmp = imwarp(A(:, :, f), tform);
            tmp = tmp(1:size(A, 1), 1:size(A, 2));
            A(:, :, f) = tmp;
            T_matrix(:, :, f-1) = tform.T;
        end
        
        x_motion = squeeze(T_matrix(3, 1, :));
        y_motion = squeeze(T_matrix(3, 2, :));
        
        
        
        % generage single averaged frame for registration
        size(x_motion)
        size(y_motion)
        
        if n == 1
            parfor f = 3 : size(data4avg, 3)
                data4avg(:, :, f) = imtranslate(data4avg(:, :, f), [x_motion(f-2)*2, y_motion(f-2)*2]);
            end
            avgFrame = mean(data4avg, 3);
            avgFrame = imrotate(avgFrame, rotateAngle);
            clear data4avg
            singleFrame_fn = [results_dir, '/','avgFrame.TIFF'];
            saveTiff16bit(avgFrame, singleFrame_fn, x_res, y_res)
        end
        
        
        % UV channel
        if n == 1
            fixed_B = B(:, :, 2);
        else
            if ~isTranslate2Initial
                fixed_B = B(:, :, 2);
            end
        end
        subMovie = B(range_r(1):range_r(2), range_c(1):range_c(2), :);
        sub_fixed = fixed_B(range_r(1):range_r(2), range_c(1):range_c(2));
        for f = 2 : size(B, 3)
            tform = imregtform(sub_fixed, subMovie(:, :, f), ...
                'translation', optimizer, metric);
            tmp = imwarp(B(:, :, f), tform);
            tmp = tmp(1:size(B, 1), 1:size(B, 2));
            B(:, :, f) = tmp;
        end
    end
    
    
    
    
    
    % rotate movies
    fprintf('Rotating movie')
    
    A = imrotate(A, rotateAngle);
    B = imrotate(B, rotateAngle);
    
    
    % reshape 3D array into space-time matrix
    sz = size(A); szZ=sz(3);
    npix = prod(sz(1:2));
    A = reshape(A, npix, szZ);
    B = reshape(B, npix, szZ);
    
    
    
    % load mask
    ROI = ReadImageJROI([data_dir, '/',mask_fn]);
    
    mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2), ...
        szall(1), szall(2));
    mask = imresize(mask, downSampleRatio, 'bilinear');
    mask = imrotate(mask, rotateAngle);
    
    [mask_r, mask_c] = find(mask > 0);
    mask_id = find(mask > 0);
    
    
    % remove slow drifting trend
    se = strel('line', hat, 0);
    
    aa = flipud(A(mask_id, 1:hat)')';
    bb = flipud(B(mask_id, 1:hat)')';
    A_sliced = [-aa + 2 * A(mask_id, 1), A(mask_id, :)];
    B_sliced = [-bb + 2 * B(mask_id, 1), B(mask_id, :)];
    
    
    fprintf('Detrending movie')
    
    parfor p = 1:length(mask_id)
        filtered1(p, :) = imtophat(A_sliced(p, :), se);
        filtered2(p, :) = imtophat(B_sliced(p, :), se);
    end
    clear A_sliced B_sliced
    A_filtered = zeros(size(A));
    B_filtered = zeros(size(B));
    
    A_filtered(mask_id, :) = filtered1(:, hat + 1 : end);
    B_filtered(mask_id, :) = filtered2(:, hat + 1 :end);
    
    clear filtered1 filtered2
    baseline_A = A - A_filtered;
    baseline_B = B - B_filtered;
    mean_A = mean(baseline_A, 2);
    mean_B = mean(baseline_B, 2);
    clear baseline_A baseline_B A B
    delete(gcp('nocreate'))
    
    
    
    % regress signal in violet frames out of signal from blue frames and
    % calculate deltaF/F.
    for w = waveLengthId
        
        clear dA
        
        switch w
            
            case 3    % regressed blue
                fprintf('Regressing uv-blue movie')
                
                small_A = A_filtered(mask_id, :) + repmat(mean_A(mask_id, :), 1, szZ);
                small_B = B_filtered(mask_id, :) + repmat(mean_B(mask_id, :), 1, szZ);
                for i = 1:length(mask_id)
                    g_pixel(i) = regress(small_A(i, :)', small_B(i, :)');
                end
                dA = zeros(size(A_filtered));
                dA(mask_id, :) = (A_filtered(mask_id, :) - repmat(g_pixel', 1, szZ) .* B_filtered(mask_id, :)) ./ ...
                    (repmat(mean_A(mask_id,:), 1, szZ));
                clear g_pixel
                
            case 1    % blue only
                fprintf('Saving blue movie')
                dA = zeros(size(A_filtered));
                dA(mask_id, :) = A_filtered(mask_id, :) ./ (repmat(mean_A(mask_id,:), 1, szZ));
                
            case 2     % UV only
                fprintf('Saving UV movie')
                dA = zeros(size(A_filtered));
                dA(mask_id, :) = B_filtered(mask_id, :) ./ (repmat(mean_B(mask_id,:), 1, szZ));
        end
        
        clear A_filtered B_filtered
        
        tmp = dA(mask_id, :);
        s0 = std(tmp(:));
        m0 = mean(tmp(:));
        movieRange = [-2*s0+m0, 5*s0+m0];
        
        save([results_dir, '/', fn(1:end-5), '_preprocessed_', waveLengthfn{w}, '.mat'], 'dA', 'movieRange','mask_id', 'sz', '-v7.3');
        
        
        % generate dF/F movie if needed
        if isStim
            generateMovie_addStim(dA, [data_dir, '/', stim_fnm], sz, movieRange, results_dir)
        end
        
    end
    toc
end
end

