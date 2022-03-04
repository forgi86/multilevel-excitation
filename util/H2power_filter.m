function h2power=h2power(Hu)
    [h,t] = impz(Hu);
    h2power = sum(h.^2);
end
