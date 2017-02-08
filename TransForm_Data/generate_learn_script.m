% 读样本
% Random shuffle 样本
% 并且将普通bin文本转成fvecs格式
% 样本格式：只有特征，没有任何描述数据大小的二进制文件
% clear;
% clc;

path = 'D:/c_sift_fv/sift_fv/x64/9/test/codefile7.bin';
learn_path = 'I:/Hashing/Code/output/fv_final/';
% n = 110463;
m = 40000;
dim = 4480;

n=size(d,2);
idx=randsample(n,m);    % 从1~n中无放回地随机选择m个数
% learn=dataset(:,idx);
learn=d(:,idx);
clear d;
learn = learn';
%% pca
pca_dim = 1024;
[pc, ~] = eigs(cov(learn), pca_dim);
fprintf('pca done!\n');
% save pca
fprintf('saving pca files...\n');
save('./fvpca_1024.mat', 'pc');
fid = fopen('./fvpca_1024.bin', 'wb');
fwrite(fid, [size(pc,2) size(pc,1)], 'int');
fwrite(fid, reshape(pc, 1, size(pc,1)*size(pc,2)),'float');
fclose(fid);

learn = learn * pc;
learn = learn';

fvecs_write(learn_path,learn);
clear learn;
