function obj=residuals_CaT(P,ident)

    P = P .*ident.sca;
    coeff=get_coeff_param(P,ident.coeff);
    
    [t,x,Ca,TR,TJ,F] = simulate_ol_scaled(coeff,ident.x0,ident.flow_profile,ident.tsim,ident.simopt);
    TRn = ident.TRn;
    Can = ident.Can;
    varC = ident.varC;
    varT = ident.varT;
    
    v_C = (eps + varC)/(varC + varT + eps);
    v_T = (eps + varT)/(varC + varT + eps);

    obj = 1/2*(((Can - Ca)'*(Can - Ca))/v_C + (TRn-TR)'*(TRn-TR)/v_T);
    
end