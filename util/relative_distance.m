function D=relative_distance(Y1,YP)
    D = zeros(size(Y1,2),1);
    for j=1:size(Y1,2)
        D(j) = norm(Y1 - YP)/norm(YP);
    end
    D = mean(D);
end
        