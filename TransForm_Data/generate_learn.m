function generate_learn(data_path,learn_path,m, d)
%% �����������˳��
% Prameters:
% In: data_path: raw data
%     m: number of data we selected
% Out: learn_path: save data after random sampling

%%
% dataset=fvecs_read(data_path);
% n=size(dataset,2);
n=size(d,2);
idx=randsample(n,m);    % ��1~n���޷Żص����ѡ��m����
% learn=dataset(:,idx);
learn=d(:,idx);
clear d;
learn = learn';
learn = double(learn);
%% pca
pca_dim = 1024;
[pc, ~] = eigs(cov(learn), pca_dim);
fprintf('pca done!\n');
% save pca
fprintf('saving pca files...\n');
save('./fvpca_c9.mat', 'pc');
fid = fopen('./fvpca_c9.bin', 'wb');
fwrite(fid, [size(pc,2) size(pc,1)], 'int');
fwrite(fid, reshape(pc, 1, size(pc,1)*size(pc,2)),'float');
fclose(fid);

learn = learn * pc;
learn = learn';

fvecs_write(learn_path,learn);
clear learn;
end
