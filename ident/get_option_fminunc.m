function o=get_option_fminunc(thetanom)
    o=optimset('LargeScale', 'off', 'TolFun', 10^-6, 'GradObj','off','TypicalX', thetanom, 'FinDiffType', 'central', 'Display', 'iter-detailed');
end
