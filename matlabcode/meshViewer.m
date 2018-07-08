clear all;close all;clc
% 1 = mesh; 2 = cross-sections; 3 = skeleton

% % CROSS-SECTIONS
% Vcs = [];
% [V, F]=reader('CSData\src_cross_sections6.txt', 2);
% Vcs = [Vcs; V];
% [V, F]=reader('CSData\src_cross_sections10.txt', 2);
% Vcs = [Vcs; V];
% [V, F]=reader('CSData\src_cross_sections14.txt', 2);
% Vcs = [Vcs; V];
% [V, F]=reader('CSData\src_cross_sections18.txt', 2);
% Vcs = [Vcs; V];

% % SKELETON
[Ss, F]=reader('src_skeleton.txt', 3);
[Sd, F]=reader('def_skeleton.txt', 3);

% % MESH
[Vs, F]=reader('src_emb_mesh.txt', 1); 
% [Vd, F]=reader('dst_gmesh.txt', 1); 


figure; hold on
% plot3(V(:, 1), V(:, 2), V(:, 3), 'b.')
% plot3(Vemb(:, 1), Vemb(:, 2), Vemb(:, 3), 'g.')
% plot3(Vcs(:, 1), Vcs(:, 2), Vcs(:, 3), 'ro')
plot3(Vs(:, 1), Vs(:, 2), Vs(:, 3), 'r.')
% plot3(Vd(:, 1), Vd(:, 2), Vd(:, 3), 'b.')
plot3(Ss(:, 1), Ss(:, 2), Ss(:, 3), 'ro')
% plot3(Sd(:, 1), Sd(:, 2), Sd(:, 3), 'bo')
axis equal
axis tight
view([0 -1 0])