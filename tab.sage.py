import sage.all
#I belive all equations from extreme rank will have to be used.

def outerFunc(N):
    for t in listOfTabToCheck:
        build=allTab[n-1][t.restrict(n-1)]
        # build.extend{map \phi to allTab[n-1][t.anti_restrict(1).rectify.standardization()]}
        # build.extend{map \phi to allTab[n-1][t.anti_restrict(1).rectify.standardization()]}
        # li=extremeRanks(t)
        allTab[n][t]=li
        grobTab[n][t]=Ideal(conjIdealOf(t)).groebner_basis()
    return allTab

def extremeRanks(t):
    n=t.size()
    li=[]
    for i2 in range(n):
#        for j1 in range(0,n+1):
        for j1 in [x-1 for x in list(dropCol(t).entries())]:
#        for j1 in [x-1 for x in flatten(list(t.conjugate())[1:])]:
            li.append(((0,-i2,j1,n),twoRank(i2,j1,t)))
        li.append(((0,-i2,n,n),twoRank(i2,n,t)))
    return li

def compare(((a1,b1,c1,d1),r1),((a2,b2,c2,d2),r2)): 
    #returns (left condition) implies (right condition); hence we want minimal elements
    return (max(a2-a1,0)+max(b2-b1,0)+max(c2-c1,0)+max(d2-d1,0)<=r2-r1) 

def fromTab(t):
    n=t.size()
    M=matrix(9, 9, [0,x12,x13,x14,x15,x16,x17,x18,x19,0,0,x23,x24,x25,x26,x27,x28,x29,0,0,0,x34,x35,x36,x37,x38,x39,0,0,0,0,x45,x46,x47,x48,x49,0,0,0,0,0,x56,x57,x58,x59,0,0,0,0,0,0,x67,x68,x69,0,0,0,0,0,0,0,x78,x79,0,0,0,0,0,0,0,0,x89,0,0,0,0,0,0,0,0,0])
    M=M[range(n),range(n)]
    for i in range(n):
        for j in range(n):
            if not(colOf(i+1,t)<colOf(j+1,t)): M[i,j]=0
    return M

def colOf(i,t): return (t.cells_containing(i)[0][1])

def actualRank(i1,i2,j1,j2,t): 
    M=fromTab(t)
    n=t.size()
    C=block_matrix([[M[range(i1,n),range(j1)],(M**2)[range(i1,n),range(j2)]],[matrix(n-i2,j1),M[range(i2,n),range(j2)]]])
    return (rank(C))

def verifyT(t):
    n=t.size()
    M=fromTab(t)
    C=block_matrix([[M,M**2],[matrix(n,n),M]])
    for ((x1,x2,j1,j2),r) in outerFunc(n)[n][t]:
        i1=-x1
        i2=-x2
        rows1=range(i1,n)
        rows2=range(n+i2,2*n)
        cols1=range(j1)
        cols2=range(n,n+j2)
        #list.extend([((i1,i2,j1,j2),conjRank(i1,i2,j1,j2,t))])
        if not (rank(C[rows1+rows2,cols1+cols2])==r): return (i1,i2,j1,j2,False)
    return True

def eqFromRank(((a,b,c,d),r),n):
    #-a (rep. -b) is the number of rows to remove from the top (resp. bottom) row block.
    #c (rep. d) is the number of columns to keep in the left (resp. right) column block.
#    M=Mbig[range(d),range(d)]
    M=Mbig[range(n),range(n)]
    C=block_matrix([[M[range(-a,d-1),range(c)],(M**2)[range(-a,d-1),range(d)]],[matrix(d-1+b,c),M[range(-b,d-1),range(d)]]])
    return [x for x in C.minors(r+1) if x!=0]

def conjIdealOf(t): 
    n=t.size()
    l=flatten([eqFromRank(x,n) for x in allTab[n][t]])
    M=Mbig[range(n),range(n)]
    l.extend([x for x in (M**3).minors(1) if x!=0])
    return(l)

def isHookOrTwo(t): 
    if (t.height()<3): return True
    if ((len(t[0])<3) or (len(t[1])==1)): return True
    return False

def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def prettyWrite(t):
    file=open(filename,'a')
    file.write('\n')
    for x in t:
        for y in x: file.write(' '+str(y))
        file.write('\n')

def dropCol(t): return Tableau((list(t.conjugate()))[1:]).conjugate() #drop the first column of tableau T

def removeFirstColUpto(i,t):
    if t.is_zero(): return t
    t=t.conjugate()
    t=t.to_list()
    l=len(t[0])
    for j in range(l):
        if ((t[0][j])<=i): t[0][j]=None
#    while (j<l and t[0][j]<=i):
#        t[0][j]=None
#        j=j+1
    t=SkewTableau(t)
    t=t.conjugate()
    return t

def rnkTab (t): return (dropCol(t).size())

def twoRank(i2,j1,t): 
    #i2=number of rows dropped from bottom row block 
    #j1=number of rows selected in left column block 
    r=rnkTab(t.restrict(j1)) #this term is block 11
    if (i2<j1): 
        t=t.anti_restrict(i2) #if i2>j1, then remove j1 +(K\cap i2)
    else: 
        t=t.anti_restrict(j1)
        t=removeFirstColUpto(i2,t)  
    t=t.rectify()
    r=r+rnkTab(t)
    return r
