import sage.all
from copy import deepcopy
#I belive all equations from extreme rank will have to be used.

n = 4
k = GF(101) 

varlist = ['dummy']
for i in range( 1, n + 1):
    for j in range( i + 1, n + 1):
        varlist.append('x' + str(i) + str(j))

A = PolynomialRing( k, varlist)
l = iter(A.gens())
x = {}
x[0] = next(l)
for i in range(1, n + 1):
    for j in range(i + 1, n+1):
        x[i,j] = next(l)

full_upper_matrix = matrix( [ x[i,j] if i < j else 0 for j in range( 1, n + 1)] for i in range(1, n + 1))

class RankCondition():
    def __init__( self, rows, cols, rank):
        assert len(rows) == len(cols)
        assert type(rows) == type((1,2))
        assert type(cols) == type((1,2))
        self.rows = rows
        self.cols = cols
        self.rank = rank
    
    def implies( self, other):
        assert len(self.rows) == len(other.rows)
        return sum( max( y - x, 0) for x, y in zip( self.rows + self.cols, other.rows + other.cols)) <= other.rank - self.rank 

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

def generic_matrix_of(tab):
    n = tab.size()
    M = deepcopy( full_upper_matrix)
    for i in range( 1, n + 1):
        for j in range( i + 1, n + 1):
            if colOf( i, tab) >= colOf( j, tab):
                M[i-1, j-1] = 0
    return M

def colOf( i, t): 
    return t.cells_containing(i)[0][1]

def rank_of_reduced_matrix( kept_rows, kept_cols, tab):
    assert len(kept_rows) == len(kept_cols)
    M = generic_matrix_of(t)
    n = t.size()
    l = []
    for i, num_rows in enumerate(kept_rows):
        if num_rows == 0: continue
        l.append([]) 
        for j, num_cols in enumerate(kept_cols):
            N = M ** ( 1 + j - i) if j >= i else matrix([ [0] * n] * n)
            l[-1].append( N[ -num_rows:, :num_cols])
    return rank( block_matrix(l))
            
def eq_from_rankCondition( condition, t):
    M = deepcopy(full_upper_matrix)
    assert len(condition.rows) < len(t[0])
    l = []
    for i, num_rows in enumerate(condition.rows):
        if num_rows == 0: continue
        l.append([])
        for j, num_cols in enumerate(condition.cols):
            N = M ** ( 1 + j - i) if j >= i else matrix([ [0] * n] * n)
            l[-1].append( N[ -num_rows:, :num_cols])

    C = block_matrix(l)
    return Ideal([x[0]]+[ x for x in C.minors( condition.rank + 1) if x != 0 ])
#    return Ideal(x for x in C.minors( condition.rank + 1) if x != 0)

def conjectured_ideal(tab):
    n = tab.size()
    k = len(t[0]) - 1
    res = []
    I = Ideal(x[0])
    for rows in product( *([[0..n]]*k)):
        for cols in product( *([[0..n]]*k)):
            print(rows,cols)
            rank = rank_of_reduced_matrix( rows, cols, tab)
            c = RankCondition( rows, cols, rank)
            if not any( x.implies(c) for x in res):
                res.append(c)
                I += Ideal(eq_from_rankCondition( c, tab))
    return Ideal(I.groebner_basis())
    

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
        if ((t[0][j])<=i): t[0][j] = None
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

