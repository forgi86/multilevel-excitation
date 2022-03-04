function [G,V,E]=get_graph(B)
    N = size(B,1);
    L = size(B,2);
    G = sparse(N,N);
    V = 1:N;
    E = [];
    for i=1:N
        for j=1:N
            if B(i,2:L) == B(j,1:L-1)
                G(i,j) = 1;
                E = [E;[i j]];
            end
        end
    end
end
