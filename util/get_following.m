function F=get_following(B)
    N = size(B,1);
    L = size(B,2);
    F = cell(N,1);
    for i=1:N
        for j=1:N
            if B(i,2:L) == B(j,1:L-1)
                F{i} = [F{i};j];
            end
        end
    end
end
               