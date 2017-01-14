% 读样本
% Random shuffle 样本
% 并且将普通bin文本转成fvecs格式
% 样本格式：只有特征，没有任何描述数据大小的二进制文件
clear;
clc;

path = 'D:/c_sift_fv/sift_fv/x64/9/test/codefile7.bin';
learn_path = 'I:/Hashing/Code/output/C9/';
n = 110463;
m = 20000;
dim = 4480;

fid = fopen(path, 'rb');
fprintf('Reading data...\n');
dataset = fread(fid, dim * n, 'float');
dataset = reshape(dataset, dim, n);
fvecs_write([learn_path 'C9_db.fvecs'], dataset);

fprintf('Random shuffle data...\n');
idx=randsample(n,m);    % 从1~n中无放回地随机选择m个数
learn=dataset(:,idx);

fvecs_write([learn_path 'C9_2w_learn.fvecs'],learn);