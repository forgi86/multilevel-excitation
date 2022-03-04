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
[P I1 opt] = experiment_design_multilevel(2,10,FM);
idx = find(P > 0.01);
P = P(idx);

T = data.T;

tt = [];
flowt = [];

t0 = 0;

tt =    [tt; t0; t0+T; t0+2*T; t0 + 3*T; t0+4*T; t0+5*T; t0+6*T; t0 + 7*T; t0+8*T];
flowt = [flowt; u1*0.6; u1*1.4; u1*1.4; u1*1.4; u1*0.6; u1*1.4; u1*1.4; u1*1.4; u1*0.6; ];

t0 = tt(end) + T;
for i=1:9
    tt = [tt; t0; t0+T; t0+2*T; t0 + 3*T;];
    flowt = [flowt; u1*1.4; u1*1.4; u1*1.4; u1*0.6];
    t0 = t0 + 4*T;
end
tt = [tt; t0+T; t0+2*T; t0+3*T; t0+4*T; t0+5*T; t0+6*T; t0+7*T; t0+8*T; t0+9*T; t0+10*T; t0+11*T; t0+12*T; t0+13*T];
flowt = [flowt; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6; u1*0.6];

flow_profile_oed = input_profile(tt, flowt);
tsim = 0:ts:tt(end);

%% THEORETICAL COVARIANCE OED %%
Ioed = 40*I1;
%% THEORETICAL COVARIANCE OED SIGNAL %%
[tm,xm,Cam,TRm,TJm,sk0m,sEm,slambdam,sUAJm,Fm] = simulate_ol_scaled_sensitivity(coeff_real,x0,flow_profile_oed,tsim,simopt);
id1 = 1;
id2 = length(tm);
sT  = [sca(1)*sk0m(id1:id2,2) sca(2)*sEm(id1:id2,2) sca(3)*slambdam(id1:id2,2) sca(4)*sUAJm(id1:id2,2)];
sC  = [sca(1)*sk0m(id1:id2,1) sca(2)*sEm(id1:id2,1) sca(3)*slambdam(id1:id2,1) sca(4)*sUAJm(id1:id2,1)];    
Ibin = sT'*sT/varT + sC'*sC/varC; % information matrix of the last experiment

%% ACTUAL COVARIANCES: MONTE CARLO SIMULATION %%
nmc=100;
optfmincon=get_option_fmincon(thetainit);
THETAOED = zeros(nmc,4);
THETARBS = THETAOED;
IRBS = cell(nmc,1);
s = RandStream('mt19937ar','Seed',3);
RandStream.setGlobalStream(s);
for j=1:nmc
%% BINARY NOISE SIMULATION %%

    [tbin,xbin,Cabin,TRbin,TJbin,Fbin] = simulate_ol_scaled(coeff_real,x0,flow_profile_oed,tsim,simopt);
    TRn = TRbin + sqrt(varT)*randn(size(TRbin));
    Can = Cabin + sqrt(varC)*randn(size(Cabin));

    ident.tsim = tsim;
    ident.flow_profile = flow_profile_oed;
    ident.varT = varT;
    ident.x0 = x0;
    ident.TRn = TRn;
    ident.Can = Can;
    ident.varC = varC;
    ident.varT = varT;
    ident.sca = sca;
    ident.coeff = coeff_real;
    ident.sca = sca;
    ident.simopt = simopt;
    resid = @(P)(residuals_CaT(P,ident));
    [thetaoed,fval,exitflag,output,lambda,grad,hessian] = fmincon(resid, thetainit,[],[],[],[],...
          [15 50 -100 5], [30 80 -60 20],[],optfmincon);
    THETAOED(j,:) = thetaoed;

    %% RBS INPUT GENERATION %%
    N = length(tt);
    flowt = idinput(N, 'rbs', [0 1],  [0.6*u1 1.4*u1]);
    flow_profile_rbs = input_profile(tt, flowt);

    %% RBS SIMULATION %%
    [t,x,Ca,TR,TJ,F] = simulate_ol_scaled(coeff_real,x0,flow_profile_rbs,tsim,simopt);
    TRn = TR + sqrt(varT)*randn(size(TR));
    Can = Ca + sqrt(varC)*randn(size(Ca));

    ident.flow_profile = flow_profile_rbs;
    ident.varT = varT;
    ident.x0 = x0;
    ident.TRn = TRn;
    ident.Can = Can;
    ident.varC = varC;
    ident.varT = varT;
    ident.sca = sca;
    ident.coeff = coeff_real;
    ident.sca = sca;
    ident.simopt = simopt;
    resid = @(P)(residuals_CaT(P,ident));
    [thetapr,fval,exitflag,output,lambda,grad,hessian] = fmincon(resid, thetainit,[],[],[],[],...
           [15 50 -100 5], [30 80 -60 20],[],optfmincon);
    THETARBS(j,:) = thetapr;
    
    %%  RBS THEORETICAL COVARIANCE%%
    simopt = simset('SrcWorkspace','current', 'Solver','ode45', 'MaxStep', ts); % for batch parallel mode
    [trbs,xrbs,Carbs,TRrbs,TJrbs,sk0rbs,sErbs,slambdarbs,sUAJrbs,Frbs] = simulate_ol_scaled_sensitivity(coeff_real,x0,flow_profile_rbs,tsim,simopt);
    id1 = 1;
    id2 = length(tm);
    sT  = [sca(1)*sk0rbs(id1:id2,2) sca(2)*sErbs(id1:id2,2) sca(3)*slambdarbs(id1:id2,2) sca(4)*sUAJrbs(id1:id2,2)];
    sC  = [sca(1)*sk0rbs(id1:id2,1) sca(2)*sErbs(id1:id2,1) sca(3)*slambdarbs(id1:id2,1) sca(4)*sUAJrbs(id1:id2,1)];    
    Irbs = sT'*sT/varT + sC'*sC/varC; % information matrix of the last experiment
    IRBS{j} = Irbs;
end

save rbs_nominal;