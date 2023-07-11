ratio_matrix = dlmread('D:\\Working issues\\Mathematical modeling\\psi.txt');
sup_ratio_matrix = dlmread('D:\\Working issues\\Mathematical modeling\\sup.txt');

figure
colormap('jet')
imagesc(ratio_matrix)
caxis([0 1])
ax1 = gca;
ax1.XTick = [];
ax1.YTick = [];
xlabel('Conversion (b), ascending')
ylabel('Fragmentation(g), ascending')
title('Prion stability')

figure
colormap('jet')
imagesc(sup_ratio_matrix)
caxis([0 1])
ax2 = gca;
ax2.XTick = [];
ax2.YTick = [];
xlabel('Conversion (b), ascending')
ylabel('Fragmentation(g), ascending')
title('Soluble Sup35 proportion')