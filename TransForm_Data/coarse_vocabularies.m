function coarse_vocabularies(yael_path,learn_data_file, raw_data_file, R_init_name, coarse_name, idx_file, K)

% Coarse quantizaton with 2 subspaces
%% Prameters:
% In: yael_path
%     file_path: data file path
%     K: number of bits of each subspace
% Out: coarse_name + '.dat': save the coarse quantization centers.
%      Rinit_path: save the coarse quantization rotetion matrix R

%%
%addpath ('~/Documents/yael/yael_v401/matlab');
learn_data = fvecs_read(learn_data_file);
learn_data = learn_data';

dim=size(learn_data,2);
niter=30;
M=2;    % number of subspaces

% tic;
%% %OPQ_NP
R_init = eye(dim);
learn_data = single(learn_data);    % convert double to float
vocabSize = 2^K;
% % add implementation of K-means

center_table_init = cell (M,1);

[centers_table_opq_np, R_opq_np , ~] = train_opq_np(learn_data,M,center_table_init,R_init,niter/2,15,vocabSize);
vocab1 = centers_table_opq_np{1}';
vocab2 = centers_table_opq_np{2}';

% save R
fvecs_write([R_init_name, '.fvecs'], R_opq_np);
save([R_init_name, '.mat'], 'R_opq_np');

% save coarse quantization centers
file = fopen([coarse_name '.dat'], 'w');
dim = size(vocab1, 1);
sz = size(vocab1, 2);
fwrite(file, dim, 'int32');
fwrite(file, sz, 'int32');
fwrite(file, vocab1, 'float');
fwrite(file, vocab2, 'float');
fclose('all');
save([coarse_name '.mat'], 'vocab1', 'vocab2');

% raw_data = fvecs_read(raw_data_file);
% raw_data = single(raw_data);

% % Cal and save idx_table (data only, no size describtion!)
% idx_table = calidx(raw_data, vocab1, vocab2);
% save([idx_file '.mat'], 'idx_table');
% idx_fout = fopen([idx_file '.dat'], 'w');
% fwrite(idx_fout, idx_table', 'int32');

learn_data = learn_data';
idx_table_part = calidx(learn_data, vocab1, vocab2);
save([idx_file '_part.mat'], 'idx_table_part');

clear learn_data;
% time=toc;
end
