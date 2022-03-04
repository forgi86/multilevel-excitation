function [M,PM]=findmarginals(G,P)

    N = length(find(G>0));
    numNodes = size(G,1);
    M = sparse(numNodes, numNodes,0);

    
    cvx_begin sdp;
    cvx_precision('low')
    %cvx_solver sdpt3
    variable PM(N,1);
    
    % marginal distributions balance
    for j=1:numNodes
        Ij = find(G(:,j) ~= 0); % indices of all the direct ancestors of the node i
        Ij = Ij';
        Pj = 0;
        for i=Ij
            k = findvariablek(G,i,j);
            Pj = Pj+PM(k)*P(i);
        end
        Pj == P(j)
    end
    
    % sum up to one    
    for i=1:numNodes
        PMi = 0;
        for j=1:numNodes
            k = findvariablek(G,i,j);
            if ~isnan(k)
                PMi = PMi + PM(k);
            end
        end
        PMi == 1
    end
    
    for i=1:N
        PM(i) >= 0
        PM(i) <= 1
    end
    
    minimize(PM'*PM)
    cvx_end
    
    for k=1:length(PM)
        [i,j] = findpositionij(G,k);
        M(i,j) = PM(k);
    end
end

            
% function M=findmarginals(Gg,P)
%     G = Gg.getmatrix;
%     [i,j,s] = find(G);
%     s(:) = NaN;
%     M = sparse(i,j,s);
%     N = length(find(G>0));
%     numNodes = size(G,1);
%     
%     cvx_begin sdp;
%     cvx_precision('medium')
%     cvx_solver sdpt3
%     variable PM(N,1);
%     for i=1:numNodes
%         if find(isnan(M(1,:)))
%             D = find(G(i,:) ~= 0); % find all direct descendant D of i
%             for j=1:length(D)
%             Dj = D(j); 
%             Aj = find(G(:,Dj) ~= 0); % all the direct ancestors of Dj
%             PDj = 0;
%             for l=1:length(A)
%                 k = findvariablek(l,j);
%                 PDj = PDj + PM(k)*P(l);
%             end
%             PDj == P(j)
%             
%             
% end
