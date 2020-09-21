clear all
close all
clc
panelA = imread('Fig7a.jpg');
panelB = imread('Fig7b.jpg');


subplot('Position',[0. 0.7 0.7 0.3])
imshow(panelA)

subplot('Position',[0.7 0.75 0.2 0.3])

text(-0.77,0.48,'(5''-UTR)','FontSize',20)

text(-0.71, 0.4,'(nsp3-F106F)','FontSize',20)

text(-0.67, 0.315,'(RdRp-P323L)','FontSize',20)
text(-0.68, 0.23,' (Spike-D614G)','FontSize',20)
text(-0.66, 0.145,'(ORF3a-Q57H/ORF3d-E14*)','FontSize',20)
axis off
subplot('Position',[0. 0 1 0.7])
imshow(panelB)
set(gcf,'PaperPosition',[0 0 16 18]);
saveas(1, 'Fig7.jpg')