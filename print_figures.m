addpath('export_fig');

load rbs_nominal;

thetaprbsavg = mean(THETAOED);

Ndata = length(t);

darkgreen = [0 0.5 0];

fontsize = 14;

%% SCATTER PLOT OF
f = figure(1);
set(f,'color','w');
h=scatter(theta0(:,1)*sca(1),theta0(:,2)*sca(2),100,'ko','filled');
set(h,'Marker','square','MarkerFaceColor','Black');

hold on;
h=plot(THETAOED(:,1)*sca(1),THETAOED(:,2)*sca(2));
set(h,'LineStyle','None', 'Marker','o', 'MarkerFaceColor', 'Green', 'MarkerEdgeColor', 'None');

h=plot(THETARBS(:,1)*sca(1),THETARBS(:,2)*sca(2));
set(h,'LineStyle','None', 'Marker','*', 'MarkerFaceColor', 'Red', 'MarkerEdgeColor', 'Red');

% h=plot(THETAOED(:,1)*sca(1),THETAOED(:,2)*sca(2));
% set(h,'LineStyle','None', 'Marker','o', 'MarkerFaceColor', 'Green', 'MarkerEdgeColor', 'None');

h=scatter(theta0(:,1)*sca(1),theta0(:,2)*sca(2),100,'ko','filled');
set(h,'Marker','square','MarkerFaceColor','Black');

l = legend('$\theta^o$', '$\hat \theta^{oed}_k$','$\hat \theta^{rbs}_k$', 'Location', 'NorthWest');
set(l,'Interpreter','Latex');
xlabel('parameter k_0','FontSize',fontsize);
ylabel('parameter E','FontSize',fontsize);
ax = gca;
set(ax,'FontSize',fontsize);
export_fig('scatter_plot', '-pdf');

Prbs = cov(THETARBS);
Poed = cov(THETAOED);
Poedth = inv(Ibin);

DETIRBS = zeros(size(IRBS));
for i=1:nmc
    DETIRBS(i) = det(IRBS{i});
end

f = figure(2);
set(f,'color','w');
plot(1:nmc, det(Ibin)*ones(size(1:nmc))*1/Ndata, 'g');
hold on;
plot(1:nmc,DETIRBS*1/Ndata, '*r');
ax = gca;
set(ax,'FontSize',fontsize);
xlabel('Monte Carlo run k','FontSize',fontsize);
l = legend('$\frac{1}{Nn}det I^{oed}$', '$\frac{1}{Nn}det I^{rbs}_k$'); 
set(l,'Interpreter','Latex','FontSize',fontsize);
export_fig('information_matrix', '-pdf');

f = figure(3);
set(f,'color','w');
plot(t,Fm/1000);
ax = gca;
set(ax,'FontSize',fontsize);
xlabel('time (s)')
ylabel('flowrate (m^3/s)','FontSize',fontsize);
legend('u^{oed}(t)','FontSize',fontsize);
export_fig('optimal_input', '-pdf');



% f = figure(4);
% set(f,'color','w');
% plot(t,Frbs/1000);
% xlabel('time (s)')
% ylabel('flowrate (m^3/s)');
% legend('u^{rbs}(t)');
% export_fig('rbs_input_1', '-pdf');


% f = figure(5);
% set(f,'color','w');
% flowt = idinput(N, 'rbs', [0 1],  [0.6*u1 1.4*u1]);
% stairs(tt,flowt/1000);
% xlabel('time (s)')
% ylabel('flowrate (m^3/s)');
% legend('u^{rbs}(t)');
% export_fig('rbs_input_2', '-pdf');
% 
% f = figure(5);
% set(f,'color','w');
% flowt = idinput(Ndata, 'rbs', [0 1],  [0.6*u1 1.4*u1]);
% stairs(tt,flowt/1000);
% xlabel('time (s)')
% ylabel('flowrate (m^3/s)');
% legend('u^{rbs}(t)');
% export_fig('rbs_input_3', '-pdf');
% 
