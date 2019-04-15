load('data/raw_resp_GH146_e51_2.mat')
spec_table = readtable('data/bnc_GH146_e51_2_LH.csv');

ROI_name = who('-file','data/raw_resp_GH146_e51_2.mat');

%  = spec_table(1); 

spec_table(1,'trl_startStk').Variables:spec_table(1,'trl_endStk').Variables
%% 
Trial_timepoints = 116;
RspTensor = zeros(length(ROI_name),size(spec_table,1));
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
                                mean(startpoint:stimstart);
    end
end

%%
%% Plot response trace to 13 different stimuli
figure(2)
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

