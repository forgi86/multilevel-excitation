close all;
clear;
s = RandStream('mt19937ar','Seed',3);
RandStream.setGlobalStream(s);

addpath('simulink');
addpath('experiment_design');
addpath('util');
addpath('ident');
addpath('export_fig');
addpath('parameters');
addpath('simulink');


%% simulation parameters for the nonlinear model %%
data = load('information_matrix_subsequences_2_10');

coeff_nominal = get_coeff_nominal;
x0 = [coeff_nominal.Cass coeff_nominal.TRss coeff_nominal.TJss]';
u1 = data.u1;
ts = 10*60; % disctere-time sampling time
simopt = simset('SrcWorkspace','current', 'Solver','ode45', 'MaxStep', ts); % for batch parallel mode

sca = [1e6 1e6 1e6 1e4];
thetanom = [coeff_nominal.k0 coeff_nominal.E coeff_nominal.lambda coeff_nominal.UAJ]./sca;
theta0 = thetanom;

coeff_real = coeff_nominal;
thetainit = [20 69.4 -69.4 8.3];

%% OPTIMAL INPUT GENERATION %%

varT = 0.1^2;
varC = 0.05^2;
FM = cell(length(data.FALL),1);
for i=1:length(FM)
    FM{i} = data.FTALL{i}/varT + data.FCALL{i}/varC;
end
[P I1 opt] = experiment_design_multilevel_trace_inv(2,10,FM);
idx = find(P > 0.01);
P = P(idx);

[G,V,E] = get_graph(data.Bin(idx,2:end));
NodeID = {};
for i=1:length(V)
    label = [char('A'+ i -1) ': ' num2str(data.Bin(idx(i),end))];
    NodeID{i} = label;
    
end
C = findcycles(biograph(G));

[Gm,M] = findmarginalsfeasibility(G,P);
GM = biograph(Gm,NodeID);
GM.ShowWeights = 'on';
view(GM);

