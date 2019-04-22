load('raw_resp_GH146_e51_2.mat')
spec_table = readtable('bnc_GH146_e51_2_LH.csv');
ROI_name = who('-file','raw_resp_GH146_e51_2.mat');

spec_table(1,'trl_startStk').Variables:spec_table(1,'trl_endStk').Variables
Trial_timepoints = 116;
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
sorted_stim_name = {'PO', 'MH04', 'MH02', 'EB04', 'EB02', 'EA04', 'EA02', ...
    'Bzald04', 'Bzald02', 'Acet04', 'Acet02', '1o3o04', '1o3o02'};
sorted_timeid_list = [];
for i = 1:length(sorted_stim_name)
    row_id = find(contains(spec_table.stim1,sorted_stim_name{i}));
    sorted_timeid_list = [sorted_timeid_list, (spec_table.trl_startStk(row_id):spec_table.trl_endStk(row_id))-1508];
end
%%
figure(4);clf;
imagesc(TraceMat(:,sorted_timeid_list));hold on;
colorbar()
for trial_j = 1:size(spec_table,1)
    startpoint = (spec_table(trial_j,'trl_startStk').Variables-1508);
    endpoint = (spec_table(trial_j,'trl_endStk').Variables-1508);
    stimstart = (spec_table(trial_j,'stim1_startStk').Variables-1508);
    stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);
    vline(startpoint - 0.5,'w')
    vline(stimstart - 0.5,'r')
    vline(stimend - 0.5,'black')
    text((stimstart + startpoint)/2, -20, sorted_stim_name{trial_j}, 'FontSize', 16)
end
for i = 1:length(z_slc_idx)
    hline(max(z_slc_idx{i})+0.5,'y:')
    text(-70, mean(z_slc_idx{i}), sprintf('slc%02d',i + 3), 'FontSize', 16)
end
axis equal tight
xticklabels([])
yticklabels([])
%%
zscore_TraceMat = zscore(TraceMat,1,2);
%% Plot of the overall zscore matrix
figure(1);clf;
imagesc(zscore_TraceMat);hold on;
colorbar()
for trial_j = 1:size(spec_table,1)
    startpoint = (spec_table(trial_j,'trl_startStk').Variables-1508);
    endpoint = (spec_table(trial_j,'trl_endStk').Variables-1508);
    stimstart = (spec_table(trial_j,'stim1_startStk').Variables-1508);
    stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);
    vline(startpoint - 0.5,'w')
    vline(stimstart - 0.5,'r')
    vline(stimend - 0.5,'black')
    text((stimstart + stimend)/2, -20, spec_table(trial_j,'stim1').Variables)
end
for i = 1:length(z_slc_idx)
    hline(max(z_slc_idx{i}))
    text(-70, mean(z_slc_idx{i}), sprintf('slc%02d',i + 3), 'FontSize', 16)
end
axis equal tight

%% Plot of the overall zscore matrix
figure(2);clf;
imagesc(zscore_TraceMat(:,sorted_timeid_list));hold on;
colorbar()
for trial_j = 1:size(spec_table,1)
    startpoint = (spec_table(trial_j,'trl_startStk').Variables-1508);
    endpoint = (spec_table(trial_j,'trl_endStk').Variables-1508);
    stimstart = (spec_table(trial_j,'stim1_startStk').Variables-1508);
    stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);
    vline(startpoint - 0.5,'w')
    vline(stimstart - 0.5,'r')
    vline(stimend - 0.5,'black')
    text((stimstart + startpoint)/2, -20, sorted_stim_name{trial_j}, 'FontSize', 16)
end
for i = 1:length(z_slc_idx)
    hline(max(z_slc_idx{i}))
    text(-70, mean(z_slc_idx{i}), sprintf('slc%02d',i + 3), 'FontSize', 16)
end
axis equal tight
xticklabels([])
yticklabels([])
%% Plot of the overall zscore matrix
Cont_subtr_zscore_TraceMat = zscore_TraceMat(:,sorted_timeid_list);
Cont_subtr_zscore_TraceMat = Cont_subtr_zscore_TraceMat - ...
    repmat(Cont_subtr_zscore_TraceMat(:,1:Trial_timepoints),1,13);
figure(3);clf;
imagesc(Cont_subtr_zscore_TraceMat);hold on;
colorbar()
for trial_j = 1:size(spec_table,1)
    startpoint = (spec_table(trial_j, 'trl_startStk').Variables-1508);
    endpoint = (spec_table(trial_j, 'trl_endStk').Variables-1508);
    stimstart = (spec_table(trial_j, 'stim1_startStk').Variables-1508);
    stimend = (spec_table(trial_j,'stim1_endStk').Variables-1508);
    vline(startpoint - 0.5,'w')
    vline(stimstart - 0.5,'r')
    vline(stimend - 0.5,'black')
    text((stimstart + stimend)/2, -20, sorted_stim_name{trial_j})
end
for i = 1:length(z_slc_idx)
    hline(max(z_slc_idx{i}))
    text(-80, mean(z_slc_idx{i}), sprintf('slice%02d',i + 3))
end
axis equal tight
%% Hierachical Clustering
odor_name = 'EB04';
row_id = find(contains(spec_table.stim1, odor_name));
startpoint = (spec_table.trl_startStk(row_id)-1508);
endpoint = (spec_table.trl_endStk(row_id)-1508);
stimstart = (spec_table.stim1_startStk(row_id)-1508);
stimend = (spec_table.stim1_endStk(row_id)-1508);
time_slice = startpoint:endpoint;
cluster_input = zscore_TraceMat(:, time_slice);

corr_mat = corrcoef(cluster_input');
% dist_mat = pdist(cluster_input');
figure(6)
imagesc(corr_mat)

Z = linkage(cluster_input, 'average', 'correlation');%'euclidean'
[~,T,perm_indx] = dendrogram(Z,0);
figure(10)
imagesc(cluster_input(perm_indx,:))
vline(stimstart + 1.5 - startpoint,'r')
vline(stimend + 1.5 - startpoint,'black')
text((stimstart + stimend)/2 - startpoint, -20, odor_name, 'FontSize', 16)
%axis equal tight
colorbar()
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
