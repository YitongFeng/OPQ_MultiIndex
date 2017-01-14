% ������
% Random shuffle ����
% ���ҽ���ͨbin�ı�ת��fvecs��ʽ
% ������ʽ��ֻ��������û���κ��������ݴ�С�Ķ������ļ�
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
idx=randsample(n,m);    % ��1~n���޷Żص����ѡ��m����
learn=dataset(:,idx);

fvecs_write([learn_path 'C9_2w_learn.fvecs'],learn);