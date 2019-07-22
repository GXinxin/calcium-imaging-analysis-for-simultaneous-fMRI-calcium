function getStimROI(data_dir, save_dir, sz)
% This function detects responding pixels to stimuli based on GLM.
%
% Required input files: 
% preprocessed dF/F, spike2 file, motion estimation from calcium imaging,
% motion estimation from simultaneous fMRI recording, estimated bleaching trend
%
% Output:
% detected ROI map, GLM beta value map, GLM results saved in mat files, and
% pixel traces within ROI saved in mat files



data_path = data_dir;
save_path = save_dir;

addpath(genpath('../code/Dependencies/'))


% initialization
isfullmodel = 1; % use full glm model
downSample = 0; 
channels = [2 13]; % channel 1: camera frame, channel 2: stim
th = 1; % threshold for spike2 stimulus onset detection
stimGap = 2; % minimum gap to recognize as two separate stimuli is 2s
sampleFreq = 1041.7; % sampling frequency of the stimulus channel on spike2 is 1041.7hz
cameraFreq = 20; % camara freq is 20hz


% load inputs
movie_fn = dir(fullfile([data_path, '/', '*_preprocessed_regress.mat'])); % preprocessed calcium movie
smr_fn = dir(fullfile([data_path, '/', '*.smr'])); % spike2 file
mri_fn = dir(fullfile([data_path, '/', 'motion_E*'])); 

fprintf('loading data...')
load([data_path, '/',movie_fn.name])
load([data_path, '/', 'trend.mat']) % estimated bleaching trend
motionMRI = load([data_path, '/', mri_fn.name]); % motion estimate from simultaneous MRI recording
load([data_path, '/', movie_fn.name(1:end-24), 'motionTrace.mat']) % motion estimate from calcium imaging
frame_num = size(dA, 2);


% load spike2 files
fHeader = fopen([data_path, '/', smr_fn.name]);
for ChannelIndex = 1:length(channels)
    [data_raw{ChannelIndex}, my_header{ChannelIndex}] = SONGetChannel(fHeader,channels(ChannelIndex),'scale');
    ChannelIndex
end


% extract stimuli onsets/offsets from spike2
isStim = (data_raw{2} > th);
stimOn = find (isStim(2:end) - isStim(1:end-1) == 1) + 1;
stimOff = find(isStim(2:end) - isStim(1:end-1) == -1) + 1;
if stimOff(1) < stimOn(1)
    stimOff = stimOff(2:end);
end

stimOn = stimOn(1:length(stimOff));

Onsets = [stimOn(1); stimOn(find(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)) + 1)];
Offsets = [stimOff(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)); stimOff(end)];

Offsets = Offsets(Offsets < sampleFreq * 630);
Onsets = Onsets(1:length(Offsets));

frameStep = data_raw{1} * sampleFreq;
frameStep = frameStep(1:2:end);



% load cortex hemisphere masks and select pixels within masks
mask_fnm = 'RoiSet.zip';
ROI = ReadImageJROI([data_path, '/', mask_fnm]);
for i = 1:length(ROI)
    roi{i} = poly2mask(ROI{i}.mnCoordinates(:, 1), ROI{i}.mnCoordinates(:, 2), sz(1), sz(2));
end

mask = roi{1} + roi{2};
if downSample
    mask = imresize(mask, .5);
    dA = imresize(reshape(dA, sz(1), sz(2), frame_num), .5);
    newSz = size(dA);
    dA = reshape(dA, newSz(1)*newSz(2), newSz(3));
end
mask_id = find(mask > 0);
img_inmask = dA(mask_id, :);


% generate box car stimlus onset trace
on_frames = [];
for i = 1 : length(Onsets)
    OnFrames = find(frameStep >= Onsets(i) - sampleFreq/cameraFreq & frameStep <= Offsets(i));
    if ~isempty(OnFrames)
        frameOnset(i) = min(OnFrames);
    end
    
    if ~isempty(find(frameStep > Offsets(i)))
        frameOffset(i) = find(frameStep > Offsets(i), 1);
    end
    
    on_frames = [on_frames, frameOnset(i) : frameOffset(i)];
end

stim = zeros(size(img_inmask, 2), 1);
stim(on_frames) = 1;



if ~isfullmodel
    glm_predictor = [stim(2:5991), trend(2:5991)', x_motion(1:5990), y_motion(1:5990)];
else
    mri_interp = interp1(10 * motionMRI(:, 1), motionMRI(:, 2:end-3), 0:6000);
    glm_predictor = [stim(2:5991), trend(2:5991)', x_motion(1:5990), y_motion(1:5990), ...
        mri_interp(1:5990, :)];
end

% run glm
fprintf('running glm...')
parfor p = 1 : length(mask_id)
    [b(p, :) dev(p) stats(p)] = glmfit(glm_predictor, img_inmask(p, 2:5991));
    pvalue(p) = stats(p).p(2);
end




if downSample
    bw_map = zeros(newSz(1), newSz(2));
else
    bw_map = zeros(sz);
end


% detect half-maximum width of beta as responding ROI
b_map = bw_map;
b_map(mask_id) = b(:, 2);

max_b = max(b_map(mask_id));
min_b = min(b_map(mask_id));

half_id = find(b_map > (max_b + min_b)/2);
bw_map(half_id) = 1;


% smooth ROI (exclude dotted regions)
cc = bwconncomp(bw_map);
bw_map_smoothed = bw_map;
for c = 1 : length(cc.PixelIdxList)
    if length(cc.PixelIdxList{c}) < 30
        bw_map_smoothed(cc.PixelIdxList{c}) = 0;
    end
end

bw_map_smoothed = imresize(bw_map_smoothed, 2, 'nearest');
bw_map_s_resp = bw_map_smoothed .* imresize(roi{1}, 2); 
bw_map_s_ctrl = bw_map_smoothed .* imresize(roi{2}, 2);


% output ROI region and beta value maps
h = figure; imagesc(bw_map_smoothed); title('ROI'); axis image;
saveas(h, [save_path, '/', movie_fn.name(1:end-24), 'glmRoi.png'])

h = figure; imagesc(b_map); title('beta values'); axis image; colorbar
saveas(h, [save_path, '/', movie_fn.name(1:end-24), 'glmBeta.png'])

% save results
save([save_path, '/', movie_fn.name(1:end-24), 'glmRoi_results.mat'], ...
    'b_map', 'bw_map', 'bw_map_smoothed', 'bw_map_s_resp', 'bw_map_s_ctrl')


% generate ROI trace
bw_map_s_resp = imresize(bw_map_s_resp, 0.5);
bw_map_s_ctrl = imresize(bw_map_s_ctrl, 0.5);
mask_id_resp = bw_map_s_resp > 0;
mask_id_ctrl = bw_map_s_ctrl > 0;
pixel_resp = dA(mask_id_resp, :);

if ~isempty(mask_id_ctrl)
    pixel_ctrl = dA(mask_id_ctrl, :);
else
    pixel_ctrl = [];
end

save([save_path, '/', movie_fn.name(1:end-24), 'pixelTrace.mat'], 'pixel_resp', 'pixel_ctrl')

