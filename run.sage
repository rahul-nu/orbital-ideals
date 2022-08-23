import datetime
from multiprocessing import Pool
execfile("tab.sage.py")

n = 6
k = GF(101) 

varlist=[]
for i in range( 1, n + 1):
    for j in range( i + 1, n + 1):
        varlist.append('x' + str(i) + str(j))

A = PolynomialRing( k, varlist)
l = iter(A.gens())
x = {}
for i in range(1, n + 1):
    for j in range(i + 1, n+1):
        x[i,j] = next(l)

M = matrix( [ x[i,j] if i < j else 0 for j in range(1, n + 1) ] for i in range(1, n + 1))

gradedListOfTab = list(StandardTableaux(0))
gradedListOfTab.extend (list(StandardTableaux(1)))

allTab=[{},{StandardTableaux(1)[0]:[x]}]
#allTab[i][t]=rank conditions on tableau t. Here i is the number of boxes in t.
for n in range(N0+1,N+1):
    gradedListOfTab.extend([[t for t in StandardTableaux(n) if len(t[0])<=3]]) 
    allTab.append({})
    for t in gradedListOfTab[n]:
        generators=allTab[n-1][t.restrict(n-1)]+flatten([eqFromRank(entry,n) for entry in extremeRanks(t)]) # + {map \phi to allTab[n-1][t.anti_restrict(1).rectify.standardization()]}
        allTab[n][t]=Ideal(generators).groebner_basis()

listOfTabToCheck=zip(listOfTabToCheck,range(1,1+len(listOfTabToCheck)))
#listOfTabToCheck=[(t,i) for (t,i) in listOfTabToCheck if not (i in done)]

def myFunc((t,i)):
    I=Ideal(conjIdealOf(t))
    J=I.groebner_basis()
    file=open(filename,'a')
    prettyWrite(t)
    file.close()

filename='allData.txt'
#file=open(filename,'a')
#file.write('dict([') 
#file.close()
for x in listOfTabToCheck:
    myFunc(x)
#pool=Pool()
#pool.map(myFunc,listOfTabToCheck)
file=open(filename,'a')
file.write('])')
file.close()
