% clear;
% clc;
% %n=3745;
% data_name='fv_c9';
% K=8;    % number of bits of each subspace in fine quantization
% M=8;    % number of subspaces of fine quantizaton
% n_choosedata = 100000; % Select how many samples to train opq
% 
% % data_path=['../../Data/' data_name];
% index_path=['I:/Hashing/Code/output/' data_name];
% mkdir(index_path);
% yael_path = 'C:/Libraries/yael438';
% vlfeat_path='C:/vlfeat-0.9.20-bin/vlfeat-0.9.20';
% 
% fid_report = fopen([index_path '/' data_name '_matlab_time.txt'], 'a');
% fprintf(fid_report, '\n%s\n', datestr(now));
tic;
generate_learn([index_path '/' data_name '_base.fvecs'],[index_path '/' data_name '_learn.fvecs'], n_choosedata, d);
time = toc;
% % fprintf(fid_report, 'Random shuffle samples time: %f \n', time);
fprintf('Random shuffle samples time: %f seconds\n', time);

tic;
coarse_vocabularies(yael_path, [index_path '/' data_name '_learn.fvecs'], [index_path '/' data_name '_base.fvecs'], ...
    [index_path '/' data_name '_rinit'],[index_path '/' data_name '_coarse'], [index_path '/' data_name '_coarse_idx'], K);
time = toc;
fprintf('Coarse quantization samples time: %f  minutes\n', time / 60);
% fprintf(fid_report, 'Coarse quantization samples time: %f \n', time);

tic;
fine_vocabularies(yael_path,vlfeat_path, [index_path '/' data_name '_learn.fvecs'], [index_path '/' data_name '_rinit.fvecs'], ...
[index_path '/' data_name '_coarse_idx'], [index_path '/' data_name '_coarse'], [index_path '/' data_name '_fine'], [index_path '/' data_name '_fine_idx'], K, M);
time = toc;
fprintf('Fine quantization samples time: %f minutes\n', time / 60);
% fprintf(fid_report, 'Fine quantization samples time: %f minutes\n', time / 60);

tic;
transform_base_query([index_path '/' data_name '_base.fvecs'], [index_path '/' data_name '_query.fvecs'], [index_path '/' data_name '_rinit.fvecs'], [index_path '/' data_name '_base_NP.fvecs'], [index_path '/' data_name '_query_NP.fvecs']);
time = toc;
fprintf('%s transform data time: %f s \n', data_name, time);
fprintf(fid_report, '%s transform data time: %f  s\n', data_name, time);


