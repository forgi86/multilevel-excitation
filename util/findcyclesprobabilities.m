function [PC]=findcyclesprobabilities(C,P)
    L = zeros(length(P), length(C));
    for i=1:size(L,1)
        for j=1:size(L,2)
            L(i,j) = length(find(C{j} == i));
        end
    end
    PC = L\P;
    PC = PC/(sum(PC));
end