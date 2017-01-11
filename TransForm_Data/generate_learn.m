function generate_learn(data_path,learn_path,m)
%% 随机打乱样本顺序
% Prameters:
% In: data_path: raw data
%     m: number of data we selected
% Out: learn_path: save data after random sampling

%%
dataset=fvecs_read(data_path);
n=size(dataset,2);
idx=randsample(n,m);    % 从1~n中无放回地随机选择m个数
learn=dataset(:,idx);

fvecs_write(learn_path,learn);
end
