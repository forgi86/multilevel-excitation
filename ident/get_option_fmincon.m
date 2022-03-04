function o=get_option_fmincon(thetanom)
    o=optimset('Algorithm', 'active-set','UseParallel','always', 'TolFun', 10^-6, 'TolX', 10^-6, 'Display', 'iter-detailed', 'GradObj','off','TypicalX', thetanom);
end
