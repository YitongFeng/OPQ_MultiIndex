function generate_learn(data_path,learn_path,m)
%% �����������˳��
% Prameters:
% In: data_path: raw data
%     m: number of data we selected
% Out: learn_path: save data after random sampling

%%
dataset=fvecs_read(data_path);
n=size(dataset,2);
idx=randsample(n,m);    % ��1~n���޷Żص����ѡ��m����
learn=dataset(:,idx);

fvecs_write(learn_path,learn);
end
