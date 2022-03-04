function [Hs,K]=scaled_filter(Hu)
    h = impulse(Hu);
    K = sum(h.^2);
    Hs = Hu/sqrt(K);   
end
