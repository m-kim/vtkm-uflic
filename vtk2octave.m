[u,v] = textread("output.vtk", "%f %f %*f", "headerlines", 5);
[x,y] = meshgrid(0:1.0/511.0:1, 0:1.0/511.0:1);
u = reshape(u, 512,512);
v = reshape(v, 512,512);
f = figure('visible', 'off')
quiver(x,y,u,v);
filename = 'screenspace.png'
print(filename)