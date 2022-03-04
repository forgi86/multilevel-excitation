function coeff = get_coeff_nominal()

    Ca0 = 8.01;
    E = 69.71*10^6;
    lambda = -69.71*10^6;
    k0 = 20.75*10^6;
    R = 8314;
    rho = 801;
    cp = 3137;
    U = 851;
    rhoJ = 1000;
    cJ = 4183;
    T0 = 294;
    TCin = 294;

    conversion = 0.95;
    Fss = 4.377e-3;    
    TRss = 350;
    
    Cass = Ca0*(1-conversion);
    kss = k0*exp(-E/((TRss*R)));
    VR = (Fss*(Ca0 - Cass)/(kss*Cass));
    D = (2*VR/pi)^(1/3);
    AJ = 2*pi*D^2;
    UAJ = U*AJ;
    Qss = (Ca0 - Cass)*Fss*(-lambda) - cp*rho*Fss*(TRss - T0);
    TJss = TRss - Qss/(U*AJ);
    VJ = 1/3*AJ;
    FJss = Qss/(cJ*(TJss-TCin)*rhoJ);
    
    
    coeff.Ca0 = Ca0;
    coeff.E = E;
    coeff.lambda = lambda;
    coeff.k0 = k0;
    coeff.R = R;
    coeff.rho = rho;
    coeff.cp = cp;
    coeff.U = U;
    coeff.rhoJ = rhoJ;
    coeff.cJ = cJ;
    coeff.T0 = T0;
    coeff.TCin = TCin;
    
    coeff.conversion = conversion;
    coeff.Fss = Fss;
    
    coeff.Cass = Cass;
    coeff.TRss = TRss;    
    coeff.TJss = TJss;
    coeff.FJss = FJss;
    
    coeff.VR = VR;
    coeff.AJ = AJ;
    coeff.UAJ = UAJ;
    coeff.kss = kss;
    coeff.TJss = TJss;
    coeff.FJss = FJss;
    
    coeff.VR = VR;
    coeff.AJ = AJ;
    coeff.VJ = VJ;
end
    