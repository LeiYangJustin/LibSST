clear all; close all; clc

V= [-30.1223 39.9038
-29.8774 -40.0965
29.8774 40.0903
30.1224 -39.9112];

E =[1 0
1 3
3 2
2 0];
E=E+1;

figure; hold on
for i = 1:size(E, 1)
plot([V(E(i, 1), 1), V(E(i, 2), 1)], [V(E(i, 1), 2), V(E(i, 2), 2)], 'b-')
end