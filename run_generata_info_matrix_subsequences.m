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
sca = [1e6 1e6 1e6 1e4];
coeff = get_coeff_nominal;
x0 = [coeff.Cass coeff.TRss coeff.TJss]';
t0 = 0;
ts = 10*60; % disctere-time sampling time
tend = 3600*20;
tsim = (0:ts:tend)';
Tb = 60*30; % 1/T = desired bandwidth
N = length(tsim);
simopt = simset('SrcWorkspace','current', 'Solver','ode45', 'MaxStep', ts); % for batch parallel mode


%% linearization around the first equilibrium %%
T1 = 350; % first operating point
Tset = T1*ones(N,1);
reference_profile = input_profile(tsim, Tset); % set point
excitation_profile = input_profile(tsim, zeros(size(tsim))); % set point
disturbance_profile = input_profile(tsim,  zeros(size(tsim))); % set point

[x1,u1,y1,dx1] = trim('testing_platform_continuous_reaction_linearization_scaled',...
    [0 T1 0]',0, [0 0 0]',... % x0, u0, y0
    2,[],[], ...              % ix, iu, iy   
    [0 0 0]', []);           % dx0, idx

T = 5*3600;
n=10;
B = PermsRep([0.6*u1 1.4*u1], n);
Bin = PermsRep([0 1], n);

FALL = cell(length(B),1);
FTALL = cell(length(B),1);
FCALL = cell(length(B),1);
varT = 0.1^2;
varC = 0.2^2;
for i=1:length(B)
    i
    b = B(i,:);
    t0 = 0;
    tsim = [];
    u = [];
    for j=1:1
        u = [u; ((b))'];
        tsim = [tsim; t0 + (0:(n-1))'*T];
        t0 = tsim(end) + T;
    end
    tsim = [tsim; t0];
    u = [u; u(end)];
    flow_profile = input_profile(tsim,u);

    simopt = simset('SrcWorkspace','current', 'Solver','ode4', 'FixedStep', ts); % for batch parallel mode    
    [t,x,Ca,TR,TJ,sk0,sE,slambda,sUAJ,F] = simulate_ol_scaled_sensitivity(coeff,x1,flow_profile,tsim,simopt);
    id1 = find(t == (n-1)*T, 1, 'first');
    id2 = find(t == n*T, 1, 'first');
      
    
    sT  = [sca(1)*sk0(id1:id2,2) sca(2)*sE(id1:id2,2) sca(3)*slambda(id1:id2,2) sca(4)*sUAJ(id1:id2,2)];
    sC  = [sca(1)*sk0(id1:id2,1) sca(2)*sE(id1:id2,1) sca(3)*slambda(id1:id2,1) sca(4)*sUAJ(id1:id2,1)];    

    
    FALL{i} = sT'*sT/varT + sC'*sC/varC;
    
    FTALL{i} = sT'*sT;
    FCALL{i} = sC'*sC;
end

save information_matrix_subsequences_2_10;
