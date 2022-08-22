import datetime
from multiprocessing import Pool
execfile("tab.sage.py")

N0=4
N=5

k=GF(101) 
#k=GF(43399) #k=GF(104729) #k=GF(3^16) #k=QQ
#x12,x13,x14,x15,x16,x17,x18,x19,x23,x24,x25,x26,x27,x28,x29,x34,x35,x36,x37,x38,x39,x45,x46,x47,x48,x49,x56,x57,x58,x59,x67,x68,x69,x78,x79,x89=k['x12,x13,x14,x15,x16,x17,x18,x19,x23,x24,x25,x26,x27,x28,x29,x34,x35,x36,x37,x38,x39,x45,x46,x47,x48,x49,x56,x57,x58,x59,x67,x68,x69,x78,x79,x89'].gens()
#x12,x13,x14,x15,x16,x17,x18,x23,x24,x25,x26,x27,x28,x34,x35,x36,x37,x38,x45,x46,x47,x48,x56,x57,x58,x67,x68,x78=k['x12,x13,x14,x15,x16,x17,x18,x23,x24,x25,x26,x27,x28,x34,x35,x36,x37,x38,x45,x46,x47,x48,x56,x57,x58,x67,x68,x78'].gens()
#Mbig=matrix(9, 9, [0,x12,x13,x14,x15,x16,x17,x18,x19,0,0,x23,x24,x25,x26,x27,x28,x29,0,0,0,x34,x35,x36,x37,x38,x39,0,0,0,0,x45,x46,x47,x48,x49,0,0,0,0,0,x56,x57,x58,x59,0,0,0,0,0,0,x67,x68,x69,0,0,0,0,0,0,0,x78,x79,0,0,0,0,0,0,0,0,x89,0,0,0,0,0,0,0,0,0])
x,x12,x23,x34,x45,x56,x67,x78,x13,x24,x35,x46,x57,x68,x14,x25,x36,x47,x58,x15,x26,x37,x48,x16,x27,x38,x17,x28,x18=k['x,x12,x23,x34,x45,x56,x67,x78,x13,x24,x35,x46,x57,x68,x14,x25,x36,x47,x58,x15,x26,x37,x48,x16,x27,x38,x17,x28,x18'].gens()
Mbig=matrix(8, 8, [0,x12,x13,x14,x15,x16,x17,x18,0,0,x23,x24,x25,x26,x27,x28,0,0,0,x34,x35,x36,x37,x38,0,0,0,0,x45,x46,x47,x48,0,0,0,0,0,x56,x57,x58,0,0,0,0,0,0,x67,x68,0,0,0,0,0,0,0,x78,0,0,0,0,0,0,0,0])

gradedListOfTab=[[],[[[1]]]]
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
