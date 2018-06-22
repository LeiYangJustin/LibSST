clear all; close all; clc

path_data = 'E:\Research\SSTsystem\MainPrj\fDataSkeleton';
addpath(path_data);
path_tool = 'E:\Research\SymmetryMaterialsCode\PatSynUII\ThirdParty\MatlabCode';
addpath(path_tool);

V = vertex_reader_with_lineskip('section_pts.txt', 0);

X = vertex_reader_with_lineskip('convex_hull.txt', 0);

Xa = vertex_reader_with_lineskip('convex_hull_after.txt', 0);

% Xa = vertex_reader_with_lineskip('polygon.txt', 0);


%%
figure; hold on
axis equal
for i =1:size(V, 1)
    plot(V(i, 1), V(i, 2), 'c.')
%     text(V(i, 1), V(i, 2), num2str(i))
end
for i =1:size(X, 1)
    plot(X(i, 1), X(i, 2), 'ro')
end
for i =1:size(X, 1)-1
    plot([X(i, 1), X(i+1, 1)], [X(i, 2), X(i+1, 2)], 'r-')
end

for i =1:size(Xa, 1)
    plot(Xa(i, 1), Xa(i, 2), 'bs', 'markerfacecolor', 'b')
end
for i =1:size(Xa, 1)-1
    plot([Xa(i, 1), Xa(i+1, 1)], [Xa(i, 2), Xa(i+1, 2)], 'b--', 'Linewidth', 1)
end
axis equal