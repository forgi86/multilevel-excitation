function v = filter_tf(H,e)
    H = tf(H);
    v =  filter(H.num{1}, H.den{1}, e);
end