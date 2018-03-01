[x,y] = meshgrid(0:0.1:2, 0:0.1:1);
for t = 125:250
  [u,v] = doublegyre(x,y,t/25.0);
  f = figure('visible','off')
  quiver(x,y,u,v);
  filename = sprintf('doublegyre-%04d.png',t);
  print(filename)
end
