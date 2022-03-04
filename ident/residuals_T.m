function obj=residuals_T(P,ident)

    P = P .*ident.sca;
    coeff=get_coeff_param(P,ident.coeff);
    
    [t,x,Ca,TR,TJ,F] = simulate_ol_scaled(coeff,ident.x0,ident.flow_profile,ident.tsim,ident.simopt);
    TRn = ident.TRn;
    Can = ident.Can;
%    varC = ident.varC;
%    varT = ident.varT;
    
%    v_C = (eps + varC)/(varC + varT + eps);
%   v_T = (eps + varT)/(varC + varT + eps);

    obj = 1/2*((TRn-TR)'*(TRn-TR));
    
end