function [t,yb,u,v,y] = simulate_cl(C,G,H,ybarsig,rsig,esig)
    simopt = simset('SrcWorkspace', 'current');
    Ts = C.Ts;
    [t,~,yb,u,v,y] = sim('platform', [esig.time(1),esig.time(end)],simopt);
end