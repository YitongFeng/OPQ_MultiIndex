% Fine quantization
% Prameters:
% In: yael_path
%     vlfeat_path
%     file_path: data file path
%     Rinit_path: coarse quantization R for init
%     coarse_name + '.dat': the coarse quantization centers.
%     K = 8: number of bits of each subspace(fixed)
%     M: number of subspaces of fine quantizaton
%     Rinit_path: save the coarse quantization rotetion matrix R
% Out: fine_name + '.dat': save the fine quantization centers.
%      
function fine_vocabularies(yael_path,vlfeat_path,file_path,Rinit_path,idx_name, coarse_name, fine_name,fine_idx_name, K,M)
%addpath ('~/Documents/yael/yael_v401/matlab');
%run('~/Downloads/software/vlfeat-0.9.20/toolbox/vl_setup');

addpath(yael_path);
run([vlfeat_path '/toolbox/vl_setup']);
all_data = fvecs_read(file_path);

% tic;
vocabSize = 2^K;
R_opq_p = fvecs_read(Rinit_path);
all_data = R_opq_p' * all_data;
% all_data = uint8(R_opq_p' * all_data);
% all_data = uint8(all_data * 10000);

load([coarse_name '.mat'], 'vocab1', 'vocab2');
load([idx_name '_part.mat'], 'idx_table_part');
idx_table_part = idx_table_part + 1;
 
i1 = idx_table_part(:, 1);
i2 = idx_table_part(:, 2);

residual = single(all_data)- single([vocab1(:,i1); vocab2(:,i2)]);
bytes_per_point = M;

% cal opq code of residuals
D = size(residual,1) / bytes_per_point; % dim of each subspace
residual_vocab = cell(bytes_per_point,1);
%fine_idx_table = zeros(size(residual, 2), M);
dist = cell(bytes_per_point,1);
niter = 30;
for m = 1:bytes_per_point
    chunk = residual(D*m-D+1:D*m,:);
        % add implementation of K-means
    [residual_vocab{m}, idx_m] = vl_kmeans(chunk, 2^K);
    %fine_idx_table(:, m) = idx_m;
    dist{m} = vl_alldist2(residual_vocab{m});          
end

% save fine quantization centers
save([fine_name '.mat'],'residual_vocab','dist');

file = fopen([fine_name '.dat'], 'w');
vocabs_count = size(residual_vocab, 1);
each_vocab_count = size(residual_vocab{1}, 2);
each_vocab_dim = size(residual_vocab{1}, 1);
fwrite(file, vocabs_count, 'int32');
fwrite(file, each_vocab_count, 'int32');
fwrite(file, each_vocab_dim, 'int32');
for i = 1:vocabs_count
    for j = 1:each_vocab_count
        a = residual_vocab{i}(:,j);
        fwrite(file, a, 'float');
    end
end

% save fine idx table
%save([fine_idx_name '.mat'], 'fine_idx_table');
% fine_idx_table = fine_idx_table - 1;
% fine_idx_fout = fopen([fine_idx_name '.dat'], 'w');
% %fwrite(fine_idx_fout, size(residual, 2), 'int32');
% %fwrite(fine_idx_fout, M, 'int32');
% fwrite(fine_idx_fout, fine_idx_table', 'int32');

% time=toc;
fclose('all');
end
