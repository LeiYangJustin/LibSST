clear all; close all; clc

path_data = 'E:\Research\SSTsystem\MainPrj\fDataSkeleton';
addpath(path_data);
path_tool = 'E:\Research\SymmetryMaterialsCode\PatSynUII\ThirdParty\MatlabCode';
addpath(path_tool);

X = vertex_reader_with_lineskip('cross_section_before.txt', 1);
X = X(:, 2:3);
F = delaunay(X);
E = make_triangle(F, X);

A = vertex_reader_with_lineskip('cross_section_after.txt', 0);


% figure; hold on
% plot(X(:, 2), X(:, 3), 'ro')
% for i = 1:size(A, 1)
%     text(X(i, 2), X(i, 3),  num2str(i))
% end
% for i = 1:size(A, 1)
%     id_incid = find(A(i, :)==1);
%     
%     i
%     id_incid
%     for j = 1:numel(id_incid)     
%         plot([X(i, 2), X(id_incid(j), 2)], [X(i, 3), X(id_incid(j), 3)], 'r-')
%     end
% end

%%
figure; hold on
axis equal
for i =1:size(X, 1)
    plot(X(i, 1), X(i, 2), 'r.')
end

% i = 112;
% text(X(i, 1), X(i, 2), num2str(i))
% i = 93;
% text(X(i, 1), X(i, 2), num2str(i))
% i = 94; 
% text(X(i, 1), X(i, 2), num2str(i))
% i = 96;
% text(X(i, 1), X(i, 2), num2str(i))


% for i = 1:size(E, 1)
%     plot([E(i, 1), E(i, 3)], [E(i, 2), E(i, 4)], 'r-.')
% end
%
for i = 1:size(A, 1)
    id_incid = find(A(i, :)==1);
    
    i
    id_incid
    for j = 1:numel(id_incid)     
        plot([X(i, 1), X(id_incid(j), 1)], [X(i, 2), X(id_incid(j), 2)], 'ko-')
    end
end