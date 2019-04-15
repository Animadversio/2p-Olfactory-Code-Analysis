xlfiledata = xlsread('bnc_GH146_e51_2_LH.csv');
load('raw_resp_GH146_e51_2.mat');
figure;
st_start = 2013;
st_end = 2028;
trial_start = 1509;
trial_end = 3016;
datasize = size(blk02slc04_roi001,2);
plot(trial_start+(1:datasize), blk02slc04_roi003);hold on;
% (st_start-trial_start)/(trial_end - trial_start)*datasize;
%plot([(st_start-trial_start)/(trial_end - trial_start)*datasize, (st_end-trial_start)/(trial_end - trial_start)*datasize],[140,140], 'linewidth',10);

A = who('blk*');
Big_Matrix = [];
for i = 1:size(A)
    temp = eval(A{i});
    Big_Matrix = [Big_Matrix;temp];
end
%Big_Matrix = cellfun(@(x) eval(x(:)), A);
[U,S,V] = svd(Big_Matrix);


