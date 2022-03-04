function [t,x,Ca,T,TJ,Tmeas,F] = simulate_cl(coeff,Fss,x0,reference_profile, excitation_profile, disturbance_profile,C,ts,tsim)
    Tss = x0(2);
    simopt = simset('SrcWorkspace','current', 'Solver','ode15', 'MaxStep', 30); % for batch parallel mode
    [t,x,Ca,T,TJ,Tmeas,F] = sim('testing_platform_continuous_reaction_cloop',tsim,simopt);
end