function P_PSD_BASIC(PSD, fr, f, ax)

hold off;
set(0, 'currentfigure', f);
set(f, 'currentaxes', ax);
hold on;

nele = size(PSD,1);
PSD = H_SMOOTHD1(PSD');
PSD = H_SMOOTHD1(PSD');
fr = linspace(fr(1), fr(end), size(PSD,1));

imagesc(fr, 1:size(PSD,1), PSD)

c1 = colorbar;
xlabel('f (Hz)');
ylabel('Lower<-Cx (Depth)->Upper');
title('PSD Difference Across Depth');

for i = 1 : nele
    labels{i} = num2str(i);
end

set(gca,'linewidth', 2, 'fontsize', 12, 'fontweight', 'bold', 'ydir', 'rev', ...
    'ylim', [1 size(PSD,1)], 'xlim', [3 200], 'ytick', linspace(1, size(PSD,1), nele), 'yticklabel', labels)

ylabel(c1, '% Diff. from Array Mean')

end