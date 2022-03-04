function [t,x,Ca,TR,TJ,F] = simulate_ol_scaled(coeff,x0,flow_profile,tsim,simopt)
    %Tss = x0(2);
    [t,x,Ca,TR,TJ,F] = sim('testing_platform_continuous_reaction_scaled',tsim,simopt);
end