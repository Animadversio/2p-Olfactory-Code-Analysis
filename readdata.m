close all; clear all;
xlfiledata = xlsread('bnc_GH146_e51_2_LH.csv');
load('raw_resp_GH146_e51_2.mat');
figure;
st_start = 2013;
st_end = 2028;
trial_start = 1509-1;
trial_end = 3016;
datasize = size(blk02slc04_roi001,2);
plot(trial_start+(1:datasize), blk02slc04_roi003);hold on;

A = who('blk*');
Big_Matrix = [];
for i = 1:size(A)
    temp = eval(A{i});
    Big_Matrix = [Big_Matrix;temp];
end
preodor = [1509:1519]-trial_start;
odor1 = [1509,1624]-trial_start;
odor2 = [1625,1740]-trial_start;

Matrix1 = Big_Matrix(:,odor1(1):odor1(2));
% subtract baseline
baseline = mean(Big_Matrix(:,preodor),2);
baseline = repmat(baseline,1,size(Matrix1,2));
Matrix1 = Matrix1 - baseline;

figure;
imagesc(Matrix1);
%Big_Matrix = cellfun(@(x) eval(x(:)), A);
[U,S,V] = svd(Matrix1);
figure;
for i = 1:4
    subplot(4,1,i);
    plotvector = i;
    plot(S(plotvector,plotvector)*V(:,plotvector));
end

