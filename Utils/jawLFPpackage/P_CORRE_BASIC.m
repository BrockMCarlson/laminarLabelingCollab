function P_CORRE_BASIC(CORRE, f, ax)

hold off;
set(0, 'currentfigure', f);
set(f, 'currentaxes', ax);
hold on;

nele = size(CORRE,1);
CORRE = H_SMOOTHD1(CORRE');
CORRE = H_SMOOTHD1(CORRE');
imagesc(CORRE)

c1 = colorbar;
xlabel('Upper<-Cx (Depth)->Lower');
ylabel('Lower<-Cx (Depth)->Upper');
title('LFP CORR Across Depth');

for i = 1 : nele
    labels{i} = num2str(i);
end

set(gca,'linewidth', 2, 'fontsize', 12, 'fontweight', 'bold', 'ydir', 'rev', ...
    'ylim', [1 size(CORRE,1)], 'xlim', [1 size(CORRE,1)], 'ytick', linspace(1, size(CORRE,1), nele), 'yticklabel', labels, ...
    'xtick', linspace(1, size(CORRE,1), nele), 'xticklabel', labels)

ylabel(c1, 'R')


end