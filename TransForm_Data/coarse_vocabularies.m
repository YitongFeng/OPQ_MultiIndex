function coarse_vocabularies(yael_path,file_path,Rinit_path, coarse_name, idx_file, K)

% Coarse quantizaton with 2 subspaces
%% Prameters:
% In: yael_path
%     file_path: data file path
%     K: number of bits of each subspace
% Out: coarse_name + '.dat': save the coarse quantization centers.
%      Rinit_path: save the coarse quantization rotetion matrix R

%%
%addpath ('~/Documents/yael/yael_v401/matlab');
all_data = fvecs_read(file_path);
all_data = all_data';

dim=size(all_data,2);
niter=30;
M=2;    % number of subspaces

tic;
%%%OPQ_NP
R_init = eye(dim);
all_data = single(all_data);    % convert double to float
vocabSize = 2^K;
% % add implementation of K-means

center_table_init = cell (M,1);

[centers_table_opq_np R_opq_np idx_table] = train_opq_np(all_data,M,center_table_init,R_init,niter/2,10,vocabSize);
vocab1 = centers_table_opq_np{1}';
vocab2 = centers_table_opq_np{2}';

% save idx_table (data only, no size describtion!)
save([idx_file '.mat'], 'idx_table');
idx_table = idx_table - 1; % let the indicies start from 0;
idx_fout = fopen([idx_file '.dat'], 'w');
fwrite(idx_fout, idx_table', 'int32');

% save R
fvecs_write(Rinit_path,R_opq_np);

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
time=toc;
end
