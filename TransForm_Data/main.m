clear;
clc;
n=7388;
data_name='fvImgDb';
K=8;    % number of bits of each subspace in fine quantization
M=8;    % number of subspaces of fine quantizaton
n_choosedata=n;

data_path=['../../data/' data_name];
index_path=['../index/' data_name];
mkdir(index_path);

yael_path = 'C:/Libraries/yael438';
vlfeat_path='C:/vlfeat-0.9.20-bin/vlfeat-0.9.20';

generate_learn([data_path '/' data_name '_base.fvecs'],[index_path '/' data_name '_learn.fvecs'],n_choosedata);

coarse_vocabularies(yael_path, [data_path '/' data_name '_base.fvecs'],[index_path '/new5/' data_name '_rinit.fvecs'],[index_path '/new5/' data_name '_coarse'], [index_path '/new5/' data_name '_coarse_idx'], 5);

fine_vocabularies(yael_path,vlfeat_path, [data_path '/' data_name '_base.fvecs'], [index_path '/new5/' data_name '_rinit.fvecs'], [index_path '/new5/' data_name '_coarse_idx'], [index_path '/new5/' data_name '_coarse'], [index_path '/new5/' data_name '_fine2'], [index_path '/new5/' data_name '_fine_idx2'], K, M);

transform_base_query([data_path '/' data_name '_base.fvecs'], [data_path '/' data_name '_query.fvecs'], [index_path '/new/' data_name '_rinit.fvecs'], [index_path '/new/' data_name '_base_NP2.fvecs'], [index_path '/new/' data_name '_query_NP2.fvecs']);



