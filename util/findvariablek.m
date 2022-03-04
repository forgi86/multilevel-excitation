function k=findvariablek(G,i,j)
    k = 0;
    for jj=1:size(G,2)
        for ii=1:size(G,1)
            if G(ii,jj) ~= 0
                k = k+1;
            end
            if i == ii && j == jj
                if G(ii,jj) == 0
                    k = NaN;
                end
                return
            end
        end
    end
    k = NaN;
    return;
end
    