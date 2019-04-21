load('raw_resp_GH146_e51_2.mat')
spec_table = readtable('bnc_GH146_e51_2_LH.csv');
ROI_name = who('-file','raw_resp_GH146_e51_2.mat');

spec_table(1,'trl_startStk').Variables:spec_table(1,'trl_endStk').Variables
%% 
z_slc_idx = {};
for z = 4:17
    idx_slc = find(contains(ROI_name, sprintf('slc%02d',z)));
    z_slc_idx{z-3} = idx_slc;
end 
%% Trace matrix 
TraceMat = zeros(length(ROI_name),length(blk02slc15_roi005));
for i = 1:length(ROI_name)
    eval(['trace = ',ROI_name{i},';']);
    TraceMat(i, :) = trace;
end
%%
zscore_TraceMat = zscore(TraceMat,1,2);
%%
figure(1);clf;
imagesc(zscore_TraceMat);hold on;
colorbar()
for trial_j = 1:size(spec_table,1)
    startpoint = (spec_table(trial_j,'trl_startStk').Variables-1508);
    endpoint = (spec_table(trial_j,'trl_endStk').Variables-1508);
    stimstart = (spec_table(trial_j,'stim1_startStk').Variables-1508);
    stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);
    vline(startpoint,'w')
    vline(stimstart,'r')
    vline(stimend,'black')
    text((stimstart + stimend)/2, -20, spec_table(trial_j,'stim1').Variables)
end
for i = 1:length(z_slc_idx)
    hline(max(z_slc_idx{i}))
end
%% 
Trial_timepoints = 116;
tic
RspTensor = zeros(length(ROI_name),size(spec_table,1));
RspTensor_OFF = zeros(length(ROI_name),size(spec_table,1));
RspTensortrace = zeros(length(ROI_name),size(spec_table,1),Trial_timepoints);
for i = 1:length(ROI_name)
    eval(['trace = ',ROI_name{i},';']);
    for trial_j = 1:size(spec_table,1)
%         trace(spec_table(trial_j,'trl_startStk').Variables:spec_table(trial_j,'trl_endStk').Variables)
        startpoint = (spec_table(trial_j,'trl_startStk').Variables-1508);
        endpoint = (spec_table(trial_j,'trl_endStk').Variables-1508);
        stimstart = (spec_table(trial_j,'stim1_startStk').Variables-1508);
        stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);

        RspTensortrace(i, trial_j, :) = trace(startpoint:endpoint);
        RspTensor(i, trial_j) = mean(trace(stimstart:stimend)) - ...
                                mean(trace(startpoint:stimstart));
        RspTensor_OFF(i, trial_j) = mean(trace(stimend:stimend+15)) - ...
                                mean(trace(startpoint:stimstart));
    end
end
toc % 34.9 s
%%
figure(2)
imagesc(RspTensor)
figure(3)
imagesc(RspTensor_OFF)
%% Plot response trace to 13 different stimuli
figure(1)
i = 10;
eval(['trace = ',ROI_name{i},';']);
for trial_j = 1:size(spec_table,1)
    subplot(13,1,trial_j)
    startpoint = (spec_table(trial_j,'trl_startStk').Variables-1508);
    endpoint = (spec_table(trial_j,'trl_endStk').Variables-1508);
    stimstart = (spec_table(trial_j,'stim1_startStk').Variables-1508);
    stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);
    plot(trace(startpoint:endpoint))
    vline(stimstart-startpoint,'b')
    vline(stimend-startpoint,'r')
    ylabel(spec_table(trial_j, 'stim1').Variables)
    if trial_j~=size(spec_table,1)
        xticklabels([]);
    end 
end
%%
save('rsptensor.mat', 'RspTensor', 'RspTensor_OFF', 'RspTensortrace')
