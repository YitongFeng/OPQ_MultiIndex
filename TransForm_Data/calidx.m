% �����ݵ���䵽�����ľ�������
% Input:
%    points vocab1 vocab2 ��������֯��eg. points��һ��λһ������������Ϊ����ά��������Ϊ������
% Output:
%    id1, id2Ϊÿ���������䵽������ľ������ĵ��±꣨��0��ʼ����
function [id1, id2] = calidx(points, vocab1, vocab2)
    disp('Cal coarse ids...');
    n_vocab = size(vocab1, 2);
    N = size(points, 2);
    
    psub1 = points(1:size(points,1) / 2, :);
    psub2 = points(size(points,1) / 2 + 1 : size(points, 1), :);
    
    vocab_norm1 = sum(vocab1.*vocab1, 1);   % get a row vector 1*256
    vocab_norm2 = sum(vocab2.*vocab2, 1);
    psub1_norm = repmat(sum(psub1 .^ 2, 1), [n_vocab, 1]) + repmat(vocab_norm1', [1, N]); % get 256*N
    psub2_norm = repmat(sum(psub2 .^ 2, 1), [n_vocab, 1]) + repmat(vocab_norm2', [1, N]); % get 256*N
    dis1 = -2.0 * vocab1' * psub1 + psub1_norm;  % get 256 * N dis
    dis2 = -2.0 * vocab2' * psub2 + psub2_norm;  % get 256 * N dis
    [~, id1] = min(dis1, [], 1);
    [~, id2] = min(dis2, [], 1);
    id1 = id1 - 1;
    id2 = id2 - 1;
    
    