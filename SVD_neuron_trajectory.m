close all; clear all;
spec_table = readtable('bnc_GH146_e51_2_LH.csv');
%xlfiledata = xlsread('bnc_GH146_e51_2_LH.csv');
load('raw_resp_GH146_e51_2.mat');
startpoint = (spec_table(:,'trl_startStk').Variables-1508);
endpoint = (spec_table(:,'trl_endStk').Variables-1508);
stimstart = (spec_table(:,'stim1_startStk').Variables-1508);
stimend = (spec_table(:,'stim1_endStk').Variables-1508);
stimname = (spec_table(:,'stim1').Variables);

figure;
st_start = 2013;
st_end = 2028;
trial_start = startpoint(1);
trial_end = endpoint(1);
datasize = size(blk02slc04_roi001,2);
plot(trial_start+(1:datasize), blk02slc04_roi003);hold on;

A = who('blk*');
Big_Matrix = [];
for i = 1:size(A)
    temp = eval(A{i});
    Big_Matrix = [Big_Matrix;temp];
end
preodor = [startpoint(1):(endpoint(1)+10)];
trial1 = [startpoint(1),endpoint(1)];
alltrial = [startpoint(1),endpoint(end)];
%odor2 = [1625,1740]-trial_start;

Matrix1 = Big_Matrix(:,alltrial(1):alltrial(2));
% subtract baseline
baseline = mean(Big_Matrix(:,preodor),2);
baseline = repmat(baseline,1,size(Matrix1,2));
Matrix1 = Matrix1 - baseline;

figure;
imagesc(Matrix1);hold on;
for i = 1:length(startpoint)
    plot([startpoint(i),startpoint(i)], get(gca, 'Ylim'),'k');
end
%Big_Matrix = cellfun(@(x) eval(x(:)), A);
[U,S,V] = svd(Matrix1);
figure;
for i = 1:4
    subplot(4,1,i);
    plotvector = i;
    plot(S(plotvector,plotvector)*V(:,plotvector));
end


% neuron trajectory

trajectoryscore = S*V';
% choose the first 3;
trajectory = trajectoryscore(1:3,:);
% cut it into trials
trialtraj = cell(length(startpoint),1);
for i = 1:length(startpoint)
    trialtraj{i} = trajectory(:,startpoint(i):endpoint(i));
end
figure;
plotstimuli = [2,5];
for i = plotstimuli%length(startpoint)
    temp_plot = trialtraj{i};
%     if i == 13
%         temp_plot = temp_plot*2;
%     end
    plot3(temp_plot(1,:),temp_plot(2,:),temp_plot(3,:)); hold on;
end
legend(stimname(plotstimuli));


% figure;
% plot3(trialtraj{1,:}',trialtraj{2,:}',trialtraj{3,:}');
% legend(stimname);


%PCA
% pca(Matrix1);






