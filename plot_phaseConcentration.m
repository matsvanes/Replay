function plot_phaseConcentration(topphase, botphase, z_top, z_bot, FOI, whichstat)
[~, maxidx] = max(mean(z_top-z_bot,1));
n=length(topphase);

cmap = flipud(brewermap(2,'RdBu'));
subplot(2,3,[1 4])
h1 = raincloud_plot(z_bot(:,maxidx)', 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
  'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
  'box_col_match', 0);
h2 = raincloud_plot(z_top(:,maxidx)', 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
  'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);


legend([h1{1} h2{1}], {'bottom percentile', 'top percentile'});
xlabel('phase concentration'); ylabel('Probability density')
box off

subplot(2,3,[2 5])
plot([1,2],[z_bot(:,maxidx), z_top(:,maxidx)],'-o','color', [.6 .6 .6]),
xlim([0.5 2.5]), hold on,
plot([1,2], [mean(z_bot(:,maxidx)), mean(z_top(:,maxidx))], '-*k', 'LineWidth',2)
xticks([1,2])
xticklabels({'bottom percentile', 'top percentile'})
ylabel('Phase concentration (a.u.)')
title('phase concentration')

subplot(2,3,3)
rho_top=[];
rho_bot=[];
% add basic vector so smallest vector is 2 and then log transform (to get
% more equal arrows)
basicvec = max([0, diff([min(z_top(:,maxidx)), 2]), diff([min(z_bot(:,maxidx)), 2])]);
logz_bot = log10(z_bot(:,maxidx) + basicvec);
logz_top = log10(z_top(:,maxidx) + basicvec); 

for k=1:n
  rho_top=[rho_top topphase{k}(maxidx,:)];
  rho_bot=[rho_bot botphase{k}(maxidx,:)];
  [u_top(1,k), v_top(1,k)] = pol2cart(circ_mean(topphase{k}(maxidx,:)'),logz_top(k));
  [u_bot(1,k), v_bot(1,k)] = pol2cart(circ_mean(botphase{k}(maxidx,:)'),logz_bot(k));
end

polarhistogram(rho_bot), hold on, polarhistogram(rho_top)
title('all phases concatenated')
subplot(2,3,6)
% set the radial limit
yl=max([max(sqrt(u_top.^2+v_top.^2)),max(sqrt(u_bot.^2+v_bot.^2))]);
x_fake=[0 yl 0 -yl];
y_fake=[yl 0 -yl 0];
h_fake=compass(x_fake,y_fake);
hold on;
% plot the vectors
h1 = compass(u_bot, v_bot); 
h2 = compass(u_top,v_top); 
h3 = compass(mean(u_bot), mean(v_bot));
h4 = compass(mean(u_top),mean(v_top)); 
% set vector properties
for k=1:n
  h1(k).Color = cmap(1,:);
  h2(k).Color = cmap(2,:);
end
h3.Color = cmap(1,:);
h3.LineWidth = 5;
h4.Color = cmap(2,:);
h4.LineWidth = 5;
set(h_fake,'Visible','off');
title('phase direction')

suptitle(sprintf('%s phase concentration difference of %s (top vs. bottom percentile)', FOI, whichstat))

