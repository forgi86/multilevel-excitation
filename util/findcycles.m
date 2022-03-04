function C=findcycles(Gg)
G = Gg.getmatrix;
C = {};
numNodes = size(G,1); 
for n = 1:numNodes
   [D,P]=graphtraverse(G,n);
   for d = D
       if G(d,n)
           C=[C;graphpred2path(P,d)];
       end
   end
   G(n,:)=0; 
end
