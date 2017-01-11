function [id1, id2] = calidx(point, vocab1, vocab2)
    n = size(vocab1, 2);
    
    psub1 = point(1:size(point,1) / 2, :);
    psub2 = point(size(point,1) / 2 + 1 : size(point, 1), :);
    
    id1 = 1;
    mindis = norm(psub1 - vocab1(:,1));
    for i = 1:n
        dis = norm(psub1 - vocab1(:, i));
        if(dis < mindis)
            mindis = dis;
            id1 = i;
        end
    end
    
    id2 = 1;
    mindis = norm(psub2 - vocab2(:,1));
    for i = 1:n
        dis = norm(psub2 - vocab2(:, i));
        if(dis < mindis)
            mindis = dis;
            id2 = i;
        end
    end
end