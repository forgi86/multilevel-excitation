function [P I opt] = experiment_design_multilevel(m,n,F)

 addpath('c:\cvx_v2');
 addpath('c:\cvx_v2\structures');
 addpath('c:\cvx_v2\lib');
 addpath('c:\cvx_v2\functions');
 addpath('c:\cvx_v2\commands');
 addpath('c:\cvx_v2\builtins');
 
 scaI = 10^-5;
 %scaI = 1;

%% LMI CONSTANTS

if length(F) ~= m^n
    error('length(F) wrong');
end
p = length(F{1});
cvx_begin sdp;
    cvx_precision('medium')
    cvx_solver sdpt3
    variable P(m^n,1);
    I = zeros(p,p);
    for i=1:length(F)
        P(i) >= 0; % 0.1/(m^n);
        I = I + F{i}*P(i); % 10^-5 scaling
    end
    %% probability distribution %%
    sum(P) == 1;
    %% stationariety %%
    for i=0:(m^(n-1)-1)
        C1 = 0;
        C2 = 0;
        for j=0:(m-1)
            C1 = C1 + P(i + j*m^(n-1) + 1);
            C2 = C2 + P(m*i + 1 + j);
        end
        C1 == C2;
    end
    I = I*scaI;
    maximize(log_det(I))
cvx_end
opt = cvx_optval;
I = I/scaI;
end