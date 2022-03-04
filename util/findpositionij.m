function [i,j]=findpositionij(G,k)
    kk = 0;
    for j=1:size(G,2)
        for i=1:size(G,1)
            if G(i,j) ~= 0
                kk = kk+1;
            end
            if kk == k
                return
            end
        end
    end
    [i,j] = NaN;
end
    