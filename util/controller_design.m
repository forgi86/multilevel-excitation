function C=controller_design(G,H,Cold,gamma)
    W1 = H;
    W2 = sqrt(gamma)*H;
    W3 = [];
    Pidid = augw(G, W1, W2, W3);
    Pidid = balreal(Pidid);
    try
        %C = hifoo(Pidid, 3, 't');
        [C,~,~]=h2syn(Pidid);
         %C = balreal(C);
         %C = balred(C,2,'truncate');
    catch err
        C = Cold;
    end
end