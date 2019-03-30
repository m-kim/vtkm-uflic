fid = fopen('build/debug/output.vtk')
num = textscan(fid, '%*s %d %*s', 1, 'headerlines', 4)
mat = textscan(fid, '%f %f %*f', num{1});
u = mat{1};
v = mat{2};
fclose(fid)

spacing = 8
[x,y] = meshgrid(0:spacing:511, 0:spacing:511);
uprime = reshape(u, 512,512);
vprime = reshape(v, 512,512);
uprime = uprime(1:spacing:end, 1:spacing:end);
vprime = vprime(1:spacing:end, 1:spacing:end);
f = figure('visible', 'off')
quiver(x,y,uprime,vprime);
filename = 'screenspace.png'
print(filename)