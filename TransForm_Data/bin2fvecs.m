clear;
clc;

path = 'D:/c_sift_fv/sift_fv/x64/9/test/codefile7.bin';
dst_path = 'I:/Hashing/Code/output/C9/codefile7.fvecs';
n = 110463;
dim = 4480;

fid = fopen(path, 'rb');
dataset = fread(fid, dim * n, 'float');
dataset = reshape(dataset, dim, n);
