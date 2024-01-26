from isogeny33 import *
from param_load import load


# ###########################################################################
#                  File containing the KummerHash class 
#     Uses our (3,3)-isogeny formulae to construct a hash function 
# ###########################################################################


class KummerHash:
    """Class for Hash function using (3,3)-isogeny"""

    def __init__(self, secp, optimal = True):
        # Initialise the parameters for security parameters secp. 
        # By default it will assume you are using optimal strategies.
        
        params = load(secp, optimal)
        self.e = params['e']
        self.p = params['p']
        self.l = params['l']
        self.Fp2 = params['Fp2']
        self.K = params['K']
        self.gens = params['gens']
        if optimal:
            self.strategy = params['strategy']
        

    def hash(self, msg):
        """ Compute the hash function using (3,3)-isogenies WITHOUT optimal strategies

            INPUT: message msg = [m1,m2,m3] where the mi = 2^l+Random(1,3^(e-1))
            OUTPUT: a tuple [a,b,c,d], which are the fundamental theta constants of the image Kummer Surface
            
        """

        TC = TripleConstantsFromThetas(self.K[1])
        
        # Construct kernel generators Q1, Q2 of the (3^e, 3^e)-isogeny from message using Four Point Ladder
        D1 = [self.gens[0], self.gens[1], self.gens[2]]
        Q1 = threeDAC(D1, msg[0], msg[1], self.l, self.K)
        # assert OnKummer(Q1, self.K[1])
        D2 = [self.gens[3], self.gens[4], self.gens[5]]
        Q2 = threeDAC(D2, msg[1], msg[2], self.l, self.K)
        # assert OnKummer(Q2, self.K[1])
        # print("4pt ladder ok")
        
        # Computing the (3^e, 3^e)-isogeny to obtain the image Kummer 

        for ee in range(self.e-1, 0, -1):

            # print(f"{ee = }")

            R = Q1
            S = Q2
            
            # Triple to get 3-torsion points, R and S will generate a (3,3)-isogeny
            for _ in range(ee):
                R = TPL(R, TC)
                S = TPL(S, TC)

            pushed_pts, image_thetas = Isogeny33ImageAndPoints(R,S,TC,[[Q1, Q2]])

            Q1 = pushed_pts[0][0] 
            Q2 = pushed_pts[0][1] 
            
            # assert OnKummer(Q1, image_thetas)
            # assert OnKummer(Q2, image_thetas)

            TC = TripleConstantsFromThetas(image_thetas)

        return image_thetas

    def hash_optimal(self, msg):
        """ Compute the hash function using (3,3)-isogenies WITH optimal strategies

            INPUT: message msg = [m1,m2,m3] where the mi = 2^l+Random(1,3^(e-1))
            OUTPUT: a tuple [a,b,c,d], which are the fundamental theta constants of the image Kummer Surface
            
        """
        TC = TripleConstantsFromThetas(self.K[1])

         # Construct kernel generators Q1, Q2 of the (3^e, 3^e)-isogeny from message using Four Point Ladder
        D1 = [self.gens[0], self.gens[1], self.gens[2]]
        Q1 = threeDAC(D1, msg[0], msg[1], self.l, self.K)
        # assert OnKummer(Q1, self.K[1])
        D2 = [self.gens[3], self.gens[4], self.gens[5]]
        Q2 = threeDAC(D2, msg[1], msg[2], self.l, self.K)
        # assert OnKummer(Q2, self.K[1])

        # Computing the (3^e, 3^e)-isogeny to obtain the image Kummer WITH optimal strategies

        pts = []
        ind = 0
        inds = []

        R = Q1
        S = Q2
            
        for row in range(1, self.e):
            while ind < self.e - row:
                pts.append([R,S])
                inds.append(ind)
                
                m = self.strategy[self.e-ind-row]

                for _ in range(m):
                    R = TPL(R, TC)
                    S = TPL(S, TC)
                ind += m

            pts, image_thetas = Isogeny33ImageAndPoints(R,S,TC,pts)
            TC = TripleConstantsFromThetas(image_thetas)

            if len(pts) > 0:
                assert len(inds) == len(pts)
                R = pts[-1][0]
                S = pts[-1][1]
                ind = inds[-1]
                pts = pts[:-1]
                inds = inds[:-1]


        return image_thetas
            
            
        
        