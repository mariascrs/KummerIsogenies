from kummer_arithmetic import *

# ###########################################################################
#     Functions to compute the (3,3)-isogeny (see Section 4 of the paper) 
# ###########################################################################


def Isogeny33Evaluate(P, coeffs):
    
    """
        Computes the image of a point under a (3,3)-isogeny with coefficients [u1,u2,u3,u4,u5]
        INPUT: - A tuple of four projective field elements [x,y,z,t]
        - Coefficients coeffs = [u1,u2,u3,u4,u5]

        OUTPUT: - The image of [x,y,z,t] under the (3,3)-isogeny defined by coeffs
        
    """
    
    # print("Starting eval...\n")
    
    P2 = Squaring(P)

    u4xy = fp2_mul(coeffs[4], fp2_mul(P[0],P[1]))
    u4zt = fp2_mul(coeffs[4], fp2_mul(P[2],P[3]))

    p1 = fp2_mul(coeffs[0], P2[0])
    p1 = fp2_add(p1, fp2_mul(coeffs[1], P2[1]))
    p1 = fp2_add(p1, fp2_mul(coeffs[2], P2[2]))
    p1 = fp2_add(p1, fp2_mul(coeffs[3], P2[3]))
    p1 = fp2_mul(P[0], p1)
    p1 = fp2_add(p1, fp2_mul(u4zt, P[1]))

    p2 = fp2_mul(coeffs[0], P2[1])
    p2 = fp2_add(p2, fp2_mul(coeffs[1], P2[0]))
    p2 = fp2_add(p2, fp2_mul(coeffs[2], P2[3]))
    p2 = fp2_add(p2, fp2_mul(coeffs[3], P2[2]))
    p2 = fp2_mul(P[1], p2)
    p2 = fp2_add(p2, fp2_mul(u4zt, P[0]))

    p3 = fp2_mul(coeffs[0], P2[2])
    p3 = fp2_add(p3, fp2_mul(coeffs[1], P2[3]))
    p3 = fp2_add(p3, fp2_mul(coeffs[2], P2[0]))
    p3 = fp2_add(p3, fp2_mul(coeffs[3], P2[1]))
    p3 = fp2_mul(P[2], p3)
    p3 = fp2_add(p3, fp2_mul(u4xy, P[3]))

    p4 = fp2_mul(coeffs[0], P2[3])
    p4 = fp2_add(p4, fp2_mul(coeffs[1], P2[2]))
    p4 = fp2_add(p4, fp2_mul(coeffs[2], P2[1]))
    p4 = fp2_add(p4, fp2_mul(coeffs[3], P2[0]))
    p4 = fp2_mul(P[3], p4)
    p4 = fp2_add(p4, fp2_mul(u4xy, P[2]))
    
    
    
    return [p1,p2,p3,p4]


def Isogeny33ImageAndPoints(R,S,TC,pts):
    """
        Given generators of a (3,3)-subgroup computes the corresponding isogeny, and pushes any points through the isogeny.

        INPUT: - Generators <R,S> of (3,3)-subgroup
        - Tripling constants TC needed to comput isogeny
        - Points pts, given as a tuple of four projective field elements [x,y,z,t], to be pushed through the isogeny

        OUTPUT: - The image of the points in pts under the (3,3)-isogeny
        - The thetas constants of the image Kummer surface  
    """

    a = TC[0][0]
    b = TC[0][1]
    c = TC[0][2]
    d = TC[0][3]

    XR = R[0]
    YR = R[1]
    ZR = R[2]
    TR = R[3]

    XS = S[0]
    YS = S[1]
    ZS = S[2]
    TS = S[3]
    
    sO = TC[1]
    hO_inv = TC[4]

    # Constructing the isogeny coefficients 
    # print(f"Constructing the isogeny coefficients ...\n")
    
    ab = fp2_mul(a,b)
    ac = fp2_mul(a,c)
    ad = fp2_mul(a,d)
    bc = fp2_mul(b,c)
    bd = fp2_mul(b,d)
    cd = fp2_mul(c,d)
    d1p = fp2_add(ab,cd) 
    d1m = fp2_sub(ab,cd)
    d2p = fp2_add(ac,bd) 
    d2m = fp2_sub(ac,bd)
    d3p = fp2_add(ad,bc)
    d3m = fp2_sub(ad,bc)
    
    D1 = fp2_mul(d1p, d1m)
    D2 = fp2_mul(d2p, d2m)
    D3 = fp2_mul(d3p, d3m)

    D12 = fp2_mul(D1, D2)
    D13 = fp2_mul(D1, D3)
    D23 = fp2_mul(D2, D3)

    hR = Hadamard(Squaring(R))
    hS = Hadamard(Squaring(S))
    hRO = FourWayMult(hR, hO_inv)
    hSO = FourWayMult(hS, hO_inv)

    HR = Hadamard(hRO)
    HS = Hadamard(hSO)

    RXX = HR[0]
    RYY = HR[1]
    RZZ = HR[2]
    RTT = HR[3]
    SXX = HS[0]
    SYY = HS[1]
    SZZ = HS[2]
    STT = HS[3]

    XYR = fp2_mul(XR,YR)
    XZR = fp2_mul(XR,ZR)
    XTR = fp2_mul(XR,TR)
    YZR = fp2_mul(YR,ZR)
    YTR = fp2_mul(YR,TR)
    ZTR = fp2_mul(ZR,TR)
    XYS = fp2_mul(XS,YS)
    XZS = fp2_mul(XS,ZS)
    XTS = fp2_mul(XS,TS)
    YZS = fp2_mul(YS,ZS)
    YTS = fp2_mul(YS,TS)
    ZTS = fp2_mul(ZS,TS)

    Rsig = fp2_mul(d1m, fp2_add(XYR, ZTR))
    Rdel = fp2_mul(d1p, fp2_sub(XYR, ZTR))
    Ssig = fp2_mul(d1m, fp2_add(XYS, ZTS))
    Sdel = fp2_mul(d1p, fp2_sub(XYS, ZTS))
    RXY = fp2_add(Rsig, Rdel)
    RZT = fp2_sub(Rsig, Rdel)
    SXY = fp2_add(Ssig, Sdel)
    SZT = fp2_sub(Ssig, Sdel)

    Rsig = fp2_mul(d2m, fp2_add(XZR, YTR))
    Rdel = fp2_mul(d2p, fp2_sub(XZR, YTR))
    Ssig = fp2_mul(d2m, fp2_add(XZS, YTS))
    Sdel = fp2_mul(d2p, fp2_sub(XZS, YTS))
    RXZ = fp2_add(Rsig, Rdel)
    RYT = fp2_sub(Rsig, Rdel)
    SXZ = fp2_add(Ssig, Sdel)
    SYT = fp2_sub(Ssig, Sdel)

    Rsig = fp2_mul(d3m, fp2_add(XTR, YZR))
    Rdel = fp2_mul(d3p, fp2_sub(XTR, YZR))
    Ssig = fp2_mul(d3m, fp2_add(XTS, YZS))
    Sdel = fp2_mul(d3p, fp2_sub(XTS, YZS))
    RXT = fp2_add(Rsig, Rdel)
    RYZ = fp2_sub(Rsig, Rdel)
    SXT = fp2_add(Ssig, Sdel)
    SYZ = fp2_sub(Ssig, Sdel)

    c2 = fp2_mul(D23, SXY)
    c2 = fp2_add(c2, fp2_mul(D13, SXZ))
    c2 = fp2_add(c2, fp2_mul(D12, SXT))

    d2 = fp2_mul(D23, RXY)
    d2 = fp2_add(d2, fp2_mul(D13, RXZ))
    d2 = fp2_add(d2, fp2_mul(D12, RXT))

    E1 = fp2_sub(fp2_mul(d2, SZT), fp2_mul(c2, RZT))
    E1 = fp2_mul(E1, D23)

    E2 = fp2_sub(fp2_mul(SXX,RYY), fp2_mul(RXX,SYY))

    E1R = fp2_mul(E1, RXX)
    E1S = fp2_mul(E1, SXX)
    E2c = fp2_mul(E2, c2)
    E2d = fp2_mul(E2, d2)

    u1 = fp2_mul(E1R, SXX)
    u1 = fp2_add(u1, u1)

    u2 = fp2_mul(E1S, RYY)
    u2 = fp2_add(u2, fp2_mul(E1R,SYY))
    tmp = fp2_add(fp2_mul(E2d,SZT), fp2_mul(E2c,RZT))
    tmp = fp2_mul(tmp, D23)
    u2 = fp2_add(u2, tmp)

    u3 = fp2_mul(E1S, RZZ)
    u3 = fp2_add(u3, fp2_mul(E1R,SZZ))
    tmp = fp2_add(fp2_mul(E2d,SYT), fp2_mul(E2c,RYT))
    tmp = fp2_mul(tmp, D13)
    u3 = fp2_add(u3, tmp)

    u4 = fp2_mul(E1S, RTT)
    u4 = fp2_add(u4, fp2_mul(E1R,STT))
    tmp = fp2_add(fp2_mul(E2d,SYZ), fp2_mul(E2c,RYZ))
    tmp = fp2_mul(tmp, D12)
    u4 = fp2_add(u4, tmp)

    u5 = fp2_mul(E2c, d2)
    u5 = fp2_add(u5, u5)

    coeffs = [u1,u2,u3,u4,u5]
    
    # Obtain the image theta constants
    # print(f"Obtaining image thetas...\n")
    
    u5ab = fp2_mul(u5, ab)
    u5cd = fp2_mul(u5, cd)

    a2 = sO[0]
    b2 = sO[1]
    c2 = sO[2]
    d2 = sO[3]
    
    aa = fp2_mul(u1, a2)
    aa = fp2_add(aa, fp2_mul(u2,b2))
    aa = fp2_add(aa, fp2_mul(u3,c2))
    aa = fp2_add(aa, fp2_mul(u4,d2))
    aa = fp2_mul(aa, a)
    aa = fp2_add(aa, fp2_mul(u5cd,b))
    
    bb = fp2_mul(u1, b2)
    bb = fp2_add(bb, fp2_mul(u2,a2))
    bb = fp2_add(bb, fp2_mul(u3,d2))
    bb = fp2_add(bb, fp2_mul(u4,c2))
    bb = fp2_mul(bb, b)
    bb = fp2_add(bb, fp2_mul(u5cd,a))

    cc = fp2_mul(u1, c2)
    cc = fp2_add(cc, fp2_mul(u2,d2))
    cc = fp2_add(cc, fp2_mul(u3,a2))
    cc = fp2_add(cc, fp2_mul(u4,b2))
    cc = fp2_mul(cc, c)
    cc = fp2_add(cc, fp2_mul(u5ab,d))

    dd = fp2_mul(u1, d2)
    dd = fp2_add(dd, fp2_mul(u2,c2))
    dd = fp2_add(dd, fp2_mul(u3,b2))
    dd = fp2_add(dd, fp2_mul(u4,a2))
    dd = fp2_mul(dd, d)
    dd = fp2_add(dd, fp2_mul(u5ab,c))
    
    image_thetas = [aa,bb,cc,dd]

    #  pushing points through

    pushed_pts = [[Isogeny33Evaluate(P[0], coeffs), Isogeny33Evaluate(P[1], coeffs)] for P in pts]

    return pushed_pts, image_thetas
    
