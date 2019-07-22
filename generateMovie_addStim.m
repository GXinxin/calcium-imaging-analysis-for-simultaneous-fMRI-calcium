function generateMovie_addStim(imgall, stim_fnm, sz, movieRange, results_dir)

% Function to add stimulus indication and generate dF/F movie of the
% whole movie or averaged response movie
%
% Inputs:
%     imgall: the preprocessed widefield dF/F matrix
%     stim_fnm: the spike2 file name with frame time stamps and stimulus
%     time stamps
%     sz: 3d size of the preprocessed matrix


addpath(genpath('../code/Dependencies/'))


rgbColors = jet(256);
roi = zeros(sz(1:2));
roi(1:50, end-50:end) = -5;

% import stimulus data
fHeader = fopen(stim_fnm);
channels = [2 13]; % channel 1: camera frame, channel 2: stim

th = 1;
stimGap = 2; % considered to be two separate stimuli if their time difference is larger than 1s
sampleFreq = 1041.7; % Hz
cameraFreq = 20;

for ChannelIndex = 1:length(channels)
    [data_raw{ChannelIndex}, my_header{ChannelIndex}] = SONGetChannel(fHeader,channels(ChannelIndex),'scale');
    ChannelIndex
end


% detect stimulus onsets and offsets
isStim = (data_raw{2} > th);
stimOn = find (isStim(2:end) - isStim(1:end-1) == 1) + 1;
stimOff = find(isStim(2:end) - isStim(1:end-1) == -1) + 1;
if stimOff(1) < stimOn(1)
    stimOff = stimOff(2:end);
end

stimOn = stimOn(1:length(stimOff));

Onsets = [stimOn(1); stimOn(find(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)) + 1)];
Offsets = [stimOff(stimOn(2:end) - stimOn(1:end-1) > (sampleFreq *stimGap)); stimOff(end)];

Offsets = Offsets(Offsets < sampleFreq * 600);
Onsets = Onsets(1:length(Offsets));
Offsets = Offsets(1:end);



% add stim marker and get average response movie
frameStep = data_raw{1} * sampleFreq;
frameStep = frameStep(1:2:end);
frameOn = -20;
frameOff = 180;


[stimMarker_x, stimMarker_y] = meshgrid(1 : 30, sz(2)-29 : sz(2));
stimMarker_id = sub2ind(sz(1:2), stimMarker_x, stimMarker_y);
t = 1;
for i = 1 : length(Onsets)
    OnFrames = find(frameStep >= Onsets(i) - sampleFreq/cameraFreq & frameStep <= Offsets(i));
    if isempty(OnFrames)
    else
        frameOnset(i) = min(OnFrames);
    end
    
    if isempty(find(frameStep > Offsets(i)))
    else
        frameOffset(i) = find(frameStep > Offsets(i), 1);
    end
    
    if exist('frameOnset', 'var') && exist('frameOffset', 'var')...
            && (length(frameOnset) >= i) && (length(frameOffset) >= i)
        imgall(stimMarker_id(:), frameOnset(i) : frameOffset(i)) = 10;
        if (i + 1 < length(Onsets)) && (frameOnset(i) + frameOn - 1 > 0)
            avgResponseM(:, :, t) = imgall(:, frameOnset(i) + frameOn - 1 : frameOnset(i) + frameOff - 1);
            t = t + 1;
        end
    end
end


% write movies: average response movie
moviefn = 'dFFMovie';
meanResp = mean(avgResponseM, 3);
meanResp = reshape(meanResp, sz(1), sz(2), size(meanResp, 2));
meanResp2 = mat2gray(meanResp, movieRange);
frStart = 1;
frEnd = size(meanResp, 3);
[I_resp, ~] = gray2ind(meanResp2, 256);
Iarr2avi(I_resp, frStart, frEnd, [results_dir, '/', moviefn, '_avgResp_addStim.avi']); % dF/F with stim marker


% write movies: whole movie
imgall = reshape(imgall, sz(1), sz(2), sz(3));
imgall2 = mat2gray(imgall, movieRange);
[I, ~] = gray2ind(imgall2, 256);
frStart = 10;
frEnd = sz(3);
Iarr2avi(I, frStart, frEnd, [results_dir, '/', moviefn, '_addStim.avi']); % dF/F with stim marker