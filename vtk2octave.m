fid = fopen('build/debug/output.vtk')
num = textscan(fid, '%*s %d %*s', 1, 'headerlines', 4)
mat = textscan(fid, '%f %f %*f', num{1});
u = mat{1};
v = mat{2};
fclose(fid)
[x,y] = meshgrid(0:1.0/511.0:1, 0:1.0/511.0:1);
u = reshape(u, 512,512);
v = reshape(v, 512,512);
f = figure('visible', 'off')
quiver(x,y,u,v);
filename = 'screenspace.png'
print(filename)