#!/usr/bin/python3

from subprocess import Popen, PIPE
DIRECTORY= './build/release/'
TYPENAME = ['Serial', 'CUDA', 'TBB']
HEIGHT = [256, 512, 1024, 2048, 4096, 8192]
WIDTH = [512,1024,2048,4096,8192,16384]
EXECNAME = 'UFLIC_rendering_'
F = open('output', 'w')

for t in TYPENAME:
    execname = DIRECTORY+EXECNAME + t
    F.writelines(t)
    for w,h in zip(WIDTH,HEIGHT):
        print([execname, str(w),str(h)])
        p = Popen([execname, str(w),str(h)], stdout=PIPE, stderr=PIPE)
        p.wait()
        out, err = p.communicate()
        print(out.decode("utf-8"))
        F.writelines(out.decode("utf-8"))

