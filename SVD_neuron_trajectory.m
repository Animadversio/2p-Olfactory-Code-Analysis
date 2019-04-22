close all; clear all;
rng(0);
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

neuron_depth = zeros(length(A),1);
neuronnumbydepth = zeros(14,1);

for z = 4:17
    idx_slc = find(contains(A, sprintf('slc%02d',z)));
    neuron_depth(idx_slc) = z;
    neuronnumbydepth(z-3) = numel(idx_slc);
end 

preodor = [startpoint(1):(endpoint(1)+10)];
trial1 = [startpoint(1),endpoint(1)];
alltrial = [startpoint(1),endpoint(end)];
%odor2 = [1625,1740]-trial_start;


sorted_stim_name = {'PO', 'MH04', 'MH02', 'EB04', 'EB02', 'EA04', 'EA02', ...
    'Bzald04', 'Bzald02', 'Acet04', 'Acet02', '1o3o04', '1o3o02'};
sorted_timeid_list = [];
for i = 1:length(sorted_stim_name)
    row_id = find(contains(spec_table.stim1,sorted_stim_name{i}));
    sorted_timeid_list = [sorted_timeid_list, (spec_table.trl_startStk(row_id):spec_table.trl_endStk(row_id))-1508];
end

Big_Matrix = Big_Matrix(:,sorted_timeid_list);

%Matrix1 = Big_Matrix(:,alltrial(1):alltrial(2));
% subtract baseline
Matrix1 = Big_Matrix;
baseline = mean(Big_Matrix(:,preodor),2);
baseline = repmat(baseline,1,size(Matrix1,2));
Matrix1 = Matrix1 - baseline;
Matrix1 = zscore(Matrix1,0,2);

figure;
imagesc(Matrix1);hold on;
for i = 1:length(startpoint)
    plot([startpoint(i),startpoint(i)], get(gca, 'Ylim'),'k');
    text((stimstart(i) + startpoint(i))/2, -20, sorted_stim_name{i}, 'FontSize', 16)
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
plotstimuli = [4,5];
for i = plotstimuli%length(startpoint)
    temp_plot = trialtraj{i};
%     if i == 13
%         temp_plot = temp_plot*2;
%     end
    plot3(temp_plot(1,:),temp_plot(2,:),temp_plot(3,:)); hold on;
end
legend(sorted_stim_name(plotstimuli));

% PCA and clustering (K-mean)
% substract the control trial
Trial_timepoints = 116;
Trialnum = 5;
Matrix2 = Matrix1 - repmat(Matrix1(:,1:Trial_timepoints),1,13);
Matrix2 = Matrix2(:,(1:Trial_timepoints)+(Trialnum-1)*Trial_timepoints);

Num_clusters = 4;
[coeff,score,latent] = pca(Matrix2);
pcascore = score(:,1:3);
clusterID = kmeans(score(:,1:10),Num_clusters);
figure;
subplot(1,2,1);
scatter3(pcascore(:,1), pcascore(:,2),pcascore(:,3), 10, clusterID);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

% evaluate cluster
subplot(1,2,2);
silhouette(Matrix2,clusterID);

% sort by clusterID
[~,I] = sort(clusterID);
clusterID = clusterID(I,:);
Sorted_Matrix = Matrix2(I,:);
Sorted_neurondepth = neuron_depth(I,:);
figure;
imagesc(Sorted_Matrix); hold on;
for i = 1:length(startpoint)
    plot([startpoint(i),startpoint(i)], get(gca, 'Ylim'),'k');
end

% plot the average response of each cluster
figure;
for i = 1:Num_clusters
    subplot(4,4,i);
    temp = Sorted_Matrix(clusterID == i,:);
    plot(mean(temp));hold on;
    title(['Cluster ID = ', mat2str(i)]);
    ylim([-10,10]);
%     xlim([1,1508]);
%     for jj = 1:length(startpoint)
%     plot([startpoint(jj),startpoint(jj)], get(gca, 'Ylim'),'k');
%     end
end


% plot each cluster by the depth
figure;
for i = 1:Num_clusters
    subplot(4,4,i);
    temp = Sorted_neurondepth(clusterID == i);
    N = histcounts(temp,3.5:17.5);
    bar(4:17,N./neuronnumbydepth');
    title(['Cluster ID = ', mat2str(i)]);
end


% figure;
% plot3(trialtraj{1,:}',trialtraj{2,:}',trialtraj{3,:}');
% legend(stimname);


%PCA
% pca(Matrix1);






