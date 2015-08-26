import os

NUMA="numactl -i all "
path=" ../inputs/"
graph="twitter"
app = "BFS"

for recur in range(2):
    for sort in range(2):
        if(sort == 1):
            makeflag ="RECURSIVE="+str(recur)+" SORT="+str(sort)+" COARSE=1 "
            start= 60619048
            cmd = makeflag + "make clean " + app
            print cmd
            os.system(cmd)
            os.system(NUMA+ "./" + app +  " -r " + str(start) + path +graph)

            makeflag ="RECURSIVE="+str(recur)+" SORT="+str(sort)+" COARSE=0 "
            start = 60622225
            cmd = makeflag + "make clean " + app
            print cmd
            os.system(cmd)
            os.system(NUMA+ "./" + app +  " -r " + str(start) + path +graph)

        else:
            makeflag ="RECURSIVE="+str(recur)+" SORT="+str(sort)+ " " 
            print makeflag
            start= 13
            cmd = makeflag + "make clean " + app
            print cmd
            os.system(cmd)
            os.system(NUMA+ "./" + app +  " -r " + str(start) + path +graph)
