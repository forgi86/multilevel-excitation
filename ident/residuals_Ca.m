function obj=residuals_Ca(P,ident,simopt)

    coeff=get_coeff_param(P);
    
    tend = ident.tend;
    TJ_profile_id = ident.TJ_profile_id;
    Qb_profile_id = ident.Qb_profile_id;
    X0 = ident.X0;
    yn = ident.yn;
    var_Ca = ident.var_Ca;
      
    TJ_profile = TJ_profile_id;
    Qb_profile = Qb_profile_id;

    [t,xsim,TJ,Qb,ysim,z,T,x]=sim('testing_platform_semibatch_reaction',[TJ_profile.time(1) TJ_profile.time(end)],simopt);

    obj = 1/2*((yn - ysim)'*(yn - ysim));
    
end