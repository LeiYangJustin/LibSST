function [V] = vertex_reader_with_lineskip(fname, numLineSkip)
V = [];
fido = fopen(fname, 'r');
tline = fgetl(fido);
for i = 1:numLineSkip
    tline = fgetl(fido);
end

while ischar(tline)
    nums = str2num(tline);
    V = [V; nums];
    tline = fgetl(fido);
end
fclose(fido);
