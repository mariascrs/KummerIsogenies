from sage.all import *
from fp_arith import *

# ###########################################################################
#  File containing functions to perform arithmetic on the Kummer surface 
# ###########################################################################

def Hadamard(P):
    """
    The 4-way Hadamard transform

	INPUTS: - A tuple of four projective field elements [x,y,z,t]
	OUTPUTS: - The Hadamard transform [x+y+z+t,x+y-z-t,x-y+z-t,x-y-z+t]		
    """
    
    T1 = fp2_add(P[0],P[1]) 
    T2 = fp2_add(P[2],P[3])
    T3 = fp2_sub(P[0],P[1])
    T4 = fp2_sub(P[2],P[3])

    return [fp2_add(T1,T2),fp2_sub(T1,T2),fp2_add(T3,T4),fp2_sub(T3,T4)]


def Squaring(P):
    """
    The 4-way squaring

	INPUTS: - A tuple of four projective field elements [x,y,z,t]
	OUTPUTS: - The coordinate-wise squarings [x^2,y^2,z^2,t^2]	
    """
    return [fp2_sqr(x) for x in P]


def FourWayMult(V, W):
    """
	The 4-way coordinate-wise multiplication of 2 tuples

	INPUTS: - Two tuples of four projective field elements
	OUTPUTS: - One 4-tuple of their coordinate-wise products
    """

    return [fp2_mul(V[0],W[0]),fp2_mul(V[1],W[1]),fp2_mul(V[2],W[2]),fp2_mul(V[3],W[3])]


def Invert4Constants(P):
    """
    The 4-way inversion in projective 3-space

	INPUTS: - One tuple of 4-elements in P^3 
	OUTPUTS: - A tuple of elements projectively equivalent to their inverses
    """

    pi1 = fp2_mul(P[2],P[3])
    pi2 = fp2_mul(pi1,P[0])	
    pi1 = fp2_mul(pi1,P[1])
    pi3 = fp2_mul(P[0],P[1])
    pi4 = fp2_mul(pi3,P[2])
    pi3 = fp2_mul(pi3,P[3])

    return [pi1,pi2,pi3,pi4]

def OnKummer(P, thetas):
    """
    Checks if a given point is on a given Kummer surface.

	INPUTS: - P: Point in P^3, represented as a 4-tuple
 			- K: Kummer surface
	OUTPUT: - boolean: true if P lies on K, false otherwise
    """
    a,b,c,d = thetas
    a2,b2,c2,d2 = Squaring(thetas)
    a4,b4,c4,d4 = Squaring([a2,b2,c2,d2])

    A2,B2,C2,D2 = Hadamard([a2,b2,c2,d2])
    p1 = fp2_mul(fp2_mul(fp2_mul(a,b),c),d)
    p2 = fp2_mul(fp2_mul(fp2_mul(A2,B2),C2),D2)

    M = fp2_mul(p1,p2)
    m1 = fp2_inv(fp2_sub(fp2_mul(a2,d2), fp2_mul(b2,c2)))
    m2 = fp2_inv(fp2_sub(fp2_mul(b2,d2), fp2_mul(a2,c2)))
    m3 = fp2_inv(fp2_sub(fp2_mul(c2,d2), fp2_mul(a2,b2)))
    E = fp2_mul(fp2_mul(fp2_mul(M,m1),m2),m3)

    F = fp2_mul(fp2_sub(fp2_sub(fp2_add(a4,d4),b4),c4),m1)
    G = fp2_mul(fp2_sub(fp2_sub(fp2_add(b4,d4),a4),c4),m2)
    H = fp2_mul(fp2_sub(fp2_sub(fp2_add(c4,d4),a4),b4),m3)

    x = P[0]
    y = P[1]
    z = P[2]
    t = P[3]

    x2,y2,z2,t2 = Squaring(P)
    x4,y4,z4,t4 = Squaring([x2,y2,z2,t2])

    T1 = fp2_add(fp2_add(fp2_add(x4,y4),z4),t4)
    T2 = fp2_mul(fp2_mul(fp2_mul(fp2_mul(E,x),y),z),t)
    T2 = fp2_add(T2, T2)
    T3 = fp2_mul(F,fp2_add(fp2_mul(x2,t2), fp2_mul(y2,z2))) 
    T4 = fp2_mul(G,fp2_add(fp2_mul(x2,z2), fp2_mul(y2,t2)))
    T5 = fp2_mul(H,fp2_add(fp2_mul(x2,y2), fp2_mul(z2,t2)))

    return fp2_add(T1, T2) == fp2_add(fp2_add(T3,T4), T5)
    




def TripleConstantsFromThetas(thetas):
    """
    Gives the constants needed to triple points on the Kummer surface as well as for the (3,3)-isogeny.

	INPUT: - thetas: The theta constants of the Kummer surface
	OUTPUT: - The tripling constants
    """
    SO = Squaring(thetas)
    HSO = Hadamard(SO)
    K3 = Invert4Constants(thetas)
    K4 = Invert4Constants(HSO)
    return [thetas, SO, HSO, K3, K4]


def normalise(P):
    """
    Normalizes Kummer points in P^3, mainly for equality and correctness checking. Sets the 4th coordinate as 1, unless there is a zero coordinate, in which case it's a 2-torsion point, where we set the 1st coordinate as 1.

	INPUT: - P: Point in P^3, represented as a 4-tuple
	OUTPUT: - Normalised point equivalent to P
    """
 
    if fp2_mul(fp2_mul(fp2_mul(P[0],P[1]),P[2]),P[3]) != 0:
        Pinv = fp2_inv(P[3])
        return [fp2_mul(P[0],Pinv),fp2_mul(P[1],Pinv),fp2_mul(P[2],Pinv), 1]
    else:
        Pinv = fp2_inv(P[0])
        return [1, fp_mul(P[1],Pinv),fp2_mul(P[2],Pinv),fp2_mul(P[3],Pinv)]


def DBL(P, K):
    """
    (Pseudo-)doubling of a Kummer point

	INPUTS: - P: Point in P^3, represented as a 4-tuple, lying on - K: Kummer surface
	OUTPUT: - [2]P 
    """
    P = Squaring(P)
    P = Hadamard(P)
    P = Squaring(P)
    P = FourWayMult(P, K[3])
    P = Hadamard(P)
    P = FourWayMult(P, K[2])

    return P


def DBLADD(P, Q, R, K):
    """
    (Pseudo-)addition of Kummer points

	INPUTS: - P,Q,R: Points in P^3, all represented as a 4-tuple, where R=(Q+-P) lying on  - K: Kummer surface
	OUTPUT: - (Q-+P)
    """


    P2 = Squaring(P)
    Q2 = Squaring(Q)

    P2 = Hadamard(P2)
    Q2 = Hadamard(Q2)

    S = FourWayMult(P2, K[3])
    Q2 = FourWayMult(Q2, S)

    Q2 = Hadamard(Q2)
    Q2 = FourWayMult(Q2, Invert4Constants(R))
    
    return Q
    

def TPL(P, TC):
    """
    (Pseudo-)tripling of a Kummer point

	INPUTS: - P: Point in P^3, represented as a 4-tuple, lying on...
 			- K: Kummer surface
	OUTPUT: - [3]P 
    """
    R = Squaring(P)
    R = Hadamard(R)
    Q = Squaring(R)
    Q = FourWayMult(Q, TC[4])
    Q = Hadamard(Q)
    Q = FourWayMult(Q, TC[3])
    Q = Squaring(Q)
    Q = Hadamard(Q)
    S = FourWayMult(R, TC[4])
    Q = FourWayMult(Q,S)
    Q = Hadamard(Q)
    P = Invert4Constants(P)
    Q = FourWayMult(Q,P)
    return Q 
    

# ###############################################################
#   3-dimensional additional chain 3DAC
# ###############################################################


def IND(I):
    """
    Indexing function used in threeDAC

    Input: an array I containing two arrays (of length two) of bits 
    Output: an index in [0,1,2,3]
    """
    L = [I[1][0]-I[0][0], I[1][1]-I[0][1]]
    if L == [-1, -1]:
        return 0
    elif L == [1, 1]:
        return 1
    elif L == [1, -1]:
        return 2
    elif L == [-1, 1]:
        return 3
    
    f"Something is wrong in IND(), L = {L}"
    assert False
    
def Encode(m,n,l):
    """
    The encoding algorithm used in threeDAC to compute the addition chain (a la DJB)

    INPUT: scalars m,n and length of chain l (here the scalars have bitsize l)
    OUTPUT: a bit first_add (determining first step of chain) and four (l-1)-bit strings b0,b1,b2,b3 (determing the other l-1 steps of the chain)
    """
    
    binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
    
    M = binary(m)
    N = binary(n)
    b0 = []
    b1 = []
    b2 = []
    b3 = []

    first_add = M[0]

    for i in range(1, l):
        mm = GF(2)(M[i-1]+M[i])
        nn = GF(2)(N[i-1]+N[i])
        mn = GF(2)(mm+nn)
        b0.append(mn) 
        b1.append(mm)
        b2.append(GF(2)(M[i]+N[i]))
        b3.append(first_add)
        first_add = GF(2)(mm + (mn+1)*first_add)
    return first_add, b0, b1, b2, b3

def DBLTHRICEADD(R, S, delRS, T, U, delTU, V, W, delVW, K):
    """
       (Pseudo)-doubling first input and (pseudo)-addition of three sets of points.

       INPUT: points R, S, delRS = (R+-S), T,U, delTU = (T+-U), V, W, delVW = (V+-W) lying on Kummer surface K
       OUTPUT: [2]R, (R+-S), (T+-U), (V+-W)
       
    """
    return DBL(R,K), DBLADD(R,S,delRS,K), DBLADD(T,U,delTU,K), DBLADD(V,W,delVW,K)

def threeDAC(Dset,m,n,l,K):
    """
        Computes P1 + [m]P2 + [n]P3, for points P1,P2,P3 given in input

        INPUT: a set D contain sets of points (and associated differences needed to compute output) on K, scalars m and n, length l of DAC, Kummer surface K 
        OUTPUT: Dset[2][0] + [m]*Dset[1][0] + [n]*Dset[1][1]
    """
    P = Dset[0]
    D = Dset[1]
    DT = Dset[2]
    
    Q = P
    
    first_add,b0,b1,b2,b3 = Encode(m,n,l) #compute the DJB chain
    
    I = [[1,1], [2,2], [1,1]]                  

    if first_add == 1:
        Q[2] = DBLADD(Q[2],D[1],D[0],K)
        I[2][1] += 1
    else:
        Q[2] = DBLADD(Q[2],D[0],D[1],K)
        I[2][0] += 1


    for i in range(l-2, -1, -1):
        
        L = [b0[i],b1[i],b2[i],b3[i]]
        
        if L == [0,0,0,0]:
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[1][0],2*I[1][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3]=DBLTHRICEADD(Q[1],Q[0],D[2],Q[2],Q[1],D[1],Q[1],Q[3],DT[IND(I)],K)
        elif L == [0,0,0,1]:
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[1][0],2*I[1][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3]=DBLTHRICEADD(Q[1],Q[0],D[2],Q[2],Q[1],D[0],Q[1],Q[3],DT[IND(I)],K)
        elif L == [0,0,1,0]:
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[1][0],2*I[1][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3] = DBLTHRICEADD(Q[1],Q[0],D[3],Q[2],Q[1],D[1],Q[1],Q[3],DT[IND(I)],K)
        elif L == [0,0,1,1]:
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[1][0],2*I[1][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3] = DBLTHRICEADD(Q[1],Q[0],D[3],Q[2],Q[1],D[0],Q[1],Q[3],DT[IND(I)],K)
        elif L == [0,1,0,0]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[0][0],2*I[0][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3] = DBLTHRICEADD(Q[0],Q[1],D[2],Q[2],Q[0],D[1],Q[0],Q[3],DT[IND(I)],K)
        elif L == [0,1,0,1]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[0][0],2*I[0][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3] = DBLTHRICEADD(Q[0],Q[1],D[2],Q[2],Q[0],D[0],Q[0],Q[3],DT[IND(I)],K)
        elif L == [0,1,1,0]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[0][0],2*I[0][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3] = DBLTHRICEADD(Q[0],Q[1],D[3],Q[2],Q[0],D[1],Q[0],Q[3],DT[IND(I)],K)
        elif L == [0,1,1,1]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[0][0],2*I[0][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[0],Q[2],Q[3] = DBLTHRICEADD(Q[0],Q[1],D[3],Q[2],Q[0],D[0],Q[0],Q[3],DT[IND(I)],K)
        elif L == [1,0,0,0]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[1],D[1],Q[0],Q[1],D[2],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,0,0,1]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[0],D[0],Q[0],Q[1],D[2],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,0,1,0]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[1],D[1],Q[0],Q[1],D[3],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,0,1,1]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[0],D[0],Q[0],Q[1],D[3],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,1,0,0]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[0],D[1],Q[0],Q[1],D[2],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,1,0,1]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[1],D[0],Q[0],Q[1],D[2],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,1,1,0]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[0][0]+I[2][0],I[0][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[0],D[1],Q[0],Q[1],D[3],Q[2],Q[3],DT[IND(I)],K)
        elif L == [1,1,1,1]: 
            I = [[I[0][0]+I[1][0],I[0][1]+I[1][1]], [2*I[2][0],2*I[2][1]], [I[1][0]+I[2][0],I[1][1]+I[2][1]]]
            Q[1],Q[2],Q[0],Q[3] = DBLTHRICEADD(Q[2],Q[1],D[0],Q[0],Q[1],D[3],Q[2],Q[3],DT[IND(I)],K)
        else:
            f"There is a problem as L = {L}"
            assert False

    return Q[3]
