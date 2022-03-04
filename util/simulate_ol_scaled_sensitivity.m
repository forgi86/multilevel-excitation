function [t,x,Ca,TR,TJ,sk0,sE,slambda,sKT,F] = simulate_ol_scaled_sensitivity(coeff,x0,flow_profile,tsim,simopt)
    %Tss = x0(2);
    [t,x,Ca,TR,TJ,sk0,sE,slambda,sKT,F] = sim('testing_platform_continuous_reaction_scaled_sensitivity',tsim,simopt);
end