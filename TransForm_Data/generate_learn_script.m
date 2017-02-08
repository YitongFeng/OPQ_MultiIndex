% ������
% Random shuffle ����
% ���ҽ���ͨbin�ı�ת��fvecs��ʽ
% ������ʽ��ֻ��������û���κ��������ݴ�С�Ķ������ļ�
% clear;
% clc;

path = 'D:/c_sift_fv/sift_fv/x64/9/test/codefile7.bin';
learn_path = 'I:/Hashing/Code/output/fv_final/';
% n = 110463;
m = 40000;
dim = 4480;

n=size(d,2);
idx=randsample(n,m);    % ��1~n���޷Żص����ѡ��m����
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
