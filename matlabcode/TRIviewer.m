clear all
close all
clc

A = vertex_reader('old_mesh.txt');
XA = A(:, 1:2);
FA = delaunay(XA);
EA = make_triangle(FA, XA);

B = vertex_reader('new_mesh.txt');
XB = B(:, 1:2);
FB = delaunay(XB);
EB = make_triangle(FB, XB);

C = vertex_reader('example_mesh.txt');
XC = C(:, 1:2);
FC = delaunay(XC);
EC = make_triangle(FC, XC);

N = nh_reader('neighborhood.txt');

%%
figure; hold on
axis equal
% for i =1:size(XA, 1)
%     plot(XA(i, 1), XA(i, 2), 'r.')
% end
% for i = 1:size(EA, 1)
%     plot([EA(i, 1), EA(i, 3)], [EA(i, 2), EA(i, 4)], 'r--')
% end

for i =1:size(XB, 1)
    plot(XB(i, 1), XB(i, 2), 'bo', 'markerfacecolor', 'b')
%     text(V(i, 1), V(i, 2), num2str(i-1))
end
for i = 1:size(EB, 1)
    plot([EB(i, 1), EB(i, 3)], [EB(i, 2), EB(i, 4)], 'b-', 'linewidth', 1)
end
% for i =1:size(N, 1)
%     nh = N{i, 1};
%     for j = 3:2:numel(nh)
%         plot([nh(1), nh(j)+nh(1)], [nh(2), nh(j+1)+nh(2)], 'm-')
%     end
% end
% figure; hold on
% axis equal
% example
for i =1:size(XC, 1)
    plot(XC(i, 1), XC(i, 2), 'go', 'markerfacecolor', 'g')
%     text(V(i, 1), V(i, 2), num2str(i-1))
end
for i = 1:size(EC, 1)
    plot([EC(i, 1), EC(i, 3)], [EC(i, 2), EC(i, 4)], 'g-', 'Linewidth', 1)
end
axis equal

%%
E= log_reader('log.txt');
figure; 
subplot(2, 1, 1)
plot([1:size(E, 1)]', E(:, 1), 's-');
grid on
subplot(2, 1, 2); hold on
plot([1:size(E, 1)]', E(:, 2)/100, 'o-');
plot([1:size(E, 1)]', E(:, 3), '*-');
hold off
grid on


% figure; hold on
% bins = [1:1:12];
% % cntX = hist(DegX, bins);
% cntV = hist(DegV, bins);
% b = bar([cntV]', 1)
% % b.FaceColor = 'Flat';
% % b(2).CData = [1, 0, 0];
% 
% % bar(bins, cntV, 'b')