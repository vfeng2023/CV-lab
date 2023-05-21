# open file
# split by ") , ("
from math import sqrt
def computeDist(p1,p2):
    distsq = 0
    for i in range(len(p1)):
        distsq += pow(p1[i]-p2[i],2)

    return sqrt(distsq)
with open("coordinates.txt","r") as f:
    lines = f.readlines()
    allcords = []
    for l in lines:
        l = l.strip()
        l = l[1:len(l)-1]
        l = l.split(") , (")
        coords = []
        for c in l:
            coords.append(list(map(float,c.split(","))))
        allcords.append(coords)
        # print(coords)

vertices = [(0,1),(1,5),(1,3),(0,2),(0,4),(4,5),(4,6),(7,5),(6,7),(6,2),(3,2),(3,7)]
for frame in allcords:
    # print(frame)
    for v1,v2 in vertices:
        print(computeDist(frame[v1],frame[v2]),end = " ")
    print()
