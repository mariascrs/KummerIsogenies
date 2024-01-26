
/////////////////////////////////////////////////////////////////////
// Useful functions for performing arithmetic on the Kummer surface
/////////////////////////////////////////////////////////////////////

function Hadamard(P); // 4A, 4S

    /////////////////////////////////////////////////////////////////
    // The four-way Hadamard transform 
    // Inputs: a tuple of four projective field elements [x,y,z,t]
    // Outputs: the Hadamard transform 
    //              [x+y+z+t,x+y-z-t,x-y+z-t,x-y-z+t]	
    /////////////////////////////////////////////////////////////////

    T1:=P[1]+P[2]; T2:=P[3]+P[4]; T3:=P[1]-P[2]; T4:=P[3]-P[4];

    return [T1+T2,T1-T2,T3+T4,T3-T4];

end function;

function Squaring(P); // 4S

    /////////////////////////////////////////////////////////////////
    // The four-way squaring
    // Inputs: a tuple of four projective field elements [x,y,z,t]
    // Outputs: the coordinate-wise squarings [x^2,y^2,z^2,t^2]	
    /////////////////////////////////////////////////////////////////

    return [P[1]^2,P[2]^2,P[3]^2,P[4]^2];

end function;

function FourWayMult(V,W) // 4M

    /////////////////////////////////////////////////////////////////
    // The four-way coordinate-wise multiplication of 2 tuples
    // Inputs: two tuples of four projective field elements
    // Outputs: one 4-tuple of their coordinate-wise products
    /////////////////////////////////////////////////////////////////

    return [V[1]*W[1],V[2]*W[2],V[3]*W[3],V[4]*W[4]];

end function;

function Invert4Constants(P) // 6M

    /////////////////////////////////////////////////////////////////
    // The four-way inversion in projective 3-space
    // Inputs: one tuple of 4-elements in P^3 
    // Outputs: a tuple of elements projectively equivalent to 
    //              their inverses
    /////////////////////////////////////////////////////////////////

    pi1:=P[3]*P[4];
    pi2:=pi1*P[1];
    pi1:=pi1*P[2];
    pi3:=P[1]*P[2];
    pi4:=pi3*P[3];
    pi3:=pi3*P[4];

    return [pi1,pi2,pi3,pi4];

end function;

function TriplingConstantsFromThetas(O)

    /////////////////////////////////////////////////////////////////
    // Gives the constants needed to triple points 
    //      on the Kummer surface as well as for the (3,3)-isogeny.
    // Inputs: the theta constants O of the Kummer surface
    // Outputs: the tripling constants
    /////////////////////////////////////////////////////////////////
    
    SO:=Squaring(O); // 4 S
    HSO:=Hadamard(SO); // 4a, 2s
    K3:=Invert4Constants(O); // 6M
    K4:=Invert4Constants(HSO); // 6M

    return [O,SO,HSO,K3,K4];  // 12M, 4a, 2s, 4S

end function;


function normalise(P) 

    /////////////////////////////////////////////////////////////////
    // Normalizes Kummer points in P^3, mainly for equality and
    // correctness checking. Sets the 4th coordinate as 1, unless 
    // there is a zero coordinate, in which case it's a 2-torsion 
    // point, where we set the 1st coordinate as 1.
    // 
    // Inputs: a point P in P^3, represented as a 4-tuple
    // Outputs: normalised point equivalent to P in P^3
    /////////////////////////////////////////////////////////////////

    if &*P ne 0 then
        return [P[1]/P[4],P[2]/P[4],P[3]/P[4],P[4]/P[4]];
    else 
        if P[4] ne 0 then 
            return [P[1]/P[4],P[2]/P[4],P[3]/P[4],P[4]/P[4]];
        elif P[3] ne 0 then 
            return [P[1]/P[3],P[2]/P[3],P[3]/P[3],P[4]/P[3]];
        elif P[2] ne 0 then
            return [P[1]/P[2],P[2]/P[2],P[3]/P[2],P[4]/P[2]];
        else
            return [P[1]/P[1],P[2]/P[1],P[3]/P[1],P[4]/P[1]];
        end if;
    end if;

end function;


function DoubleKummer(P,K); 

    /////////////////////////////////////////////////////////////////
    // (Pseudo-)doubling of a Kummer point
    // Inputs: a point P in P^3, represented as a 4-tuple, lying on
    //         Kummer surface K
    // Outputs: [2]P
    /////////////////////////////////////////////////////////////////

    P:=Squaring(P);
    P:=Hadamard(P);
    P:=Squaring(P);
    P:=FourWayMult(P,K[4]);
    P:=Hadamard(P);
    P:=FourWayMult(P,K[3]);

    //assert OnKummer(P,K);

    return normalise(P);

end function;

function PseudoAddKummer(P,Q,R,K);

    /////////////////////////////////////////////////////////////////
    // (Pseudo-)addition of Kummer points
    // Inputs: P,Q,R: Points in P^3, all represented as a 4-tuple, 
    //         where R=(Q+-P) lying on Kummer surface K
    // Outputs: (Q-+P)
    /////////////////////////////////////////////////////////////////

    P:=Squaring(P);
    Q:=Squaring(Q);
    P:=Hadamard(P);
    Q:=Hadamard(Q);
    S:=FourWayMult(P,K[4]);
    Q:=FourWayMult(Q,S);
    Q:=Hadamard(Q);
    Q:=FourWayMult(Q,Invert4Constants(R));
    
    //assert OnKummer(Q,K);

    return Q;

end function;

function TripleKummer(P,TC) //26M+12S+16a+16s

    /////////////////////////////////////////////////////////////////
    // (Pseudo-)tripling of a Kummer point
    // Inputs: a point P in P^3, represented as a 4-tuple 
    //         lying on Kummer surface K
    // Outputs: [3]P
    /////////////////////////////////////////////////////////////////


    R:=Squaring(P);                         //4S
    R:=Hadamard(R);                         //4a+4s              
    Q:=Squaring(R);                         //4S
    Q:=FourWayMult(Q,TC[5]);                //4M
    Q:=Hadamard(Q);                         //4a+4s    
    Q:=FourWayMult(Q,TC[4]);                //4M
    Q:=Squaring(Q);                         //4S
    Q:=Hadamard(Q);                         //4a+4s    
    S:=FourWayMult(R,TC[5]);                //4M
    Q:=FourWayMult(Q,S);                    //4M
    Q:=Hadamard(Q);                         //4a+4s    
    P:=Invert4Constants(P);                 //6M
    Q:=FourWayMult(Q,P);                    //4M

    return Q;                    

end function;


///////////////////////////////////////////////////////////////
// Functions for the 3-dimensional additional chain 3DAC
///////////////////////////////////////////////////////////////


function IND(I)

    ////////////////////////////////////////////////////////////////////
    // Indexing function used in threeDAC
    // Inputs: an array I containing two arrays (of length two) of bits 
    // Outputs: an index \in [0,1,2,3]
    ////////////////////////////////////////////////////////////////////
        
    case [I[2][1]-I[1][1], I[2][2]-I[1][2]]: 

        when [-1,-1]: ind:=1;
        when [1 , 1]: ind:=2;
        when [1, -1]: ind:=3;
        when [-1, 1]: ind:=4;

    end case;

    return ind;

end function;

function ENCODE(m,n,len)

    ////////////////////////////////////////////////////////////////////
    // The encoding algorithm used in threeDAC to compute the 
    //          addition chain (a la DJB)
    // Inputs: scalars m,n and length of chain l (here the scalars 
    //         have bitsize l)
    // Outputs: a bit first_add (determining first step of chain) and 
    //          four (l-1)-bit strings b0,b1,b2,b3 (determing the 
    //          other l-1 steps of the chain)
    ////////////////////////////////////////////////////////////////////

    M:=[GF(2)!i: i in IntegerToSequence(m,2)];  //bit sequence of m
    N:=[GF(2)!i: i in IntegerToSequence(n,2)];  //bit sequence of n

    b0:=[];  b1:=[]; b2:=[]; b3:=[];    
    
    first_add:=M[1];                
    for i:=1 to (len-1) do              
        mm:=M[i]+M[i+1];
        nn:=N[i]+N[i+1];
        mn:=mm+nn;
        b0[i]:=mn;
        b1[i]:=mm;
        b2[i]:=M[i+1]+N[i+1];
        b3[i]:=first_add;
        first_add:=mm+(mn+1)*first_add;     
    end for;                    

    return first_add,b0,b1,b2,b3;

end function;

function DoubleThriceAdd(R,S,delRS,T,U,delTU,V,W,delVW,K);

    ////////////////////////////////////////////////////////////////////
    // (Pseudo)-doubling first input and (pseudo)-addition of three 
    //      sets of points.
    // Inputs: points R, S, delRS = (R+-S), T,U, delTU = (T+-U), V, W, 
    //         delVW = (V+-W) lying on Kummer surface K
    // Outputs: [2]R, (R+-S), (T+-U), (V+-W)
    ////////////////////////////////////////////////////////////////////

    return DoubleKummer(R,K),
           PseudoAddKummer(R,S,delRS,K),
           PseudoAddKummer(T,U,delTU,K),
           PseudoAddKummer(V,W,delVW,K);

end function;

function threeDAC(Dset,a,b,len,K)

    ////////////////////////////////////////////////////////////////////
    // Computes P1 + [m]P2 + [n]P3, for points P1,P2,P3 given in input
    // Inputs: a set D contain sets of points (and associated differences
    //         needed to compute output) on K, scalars m and n, 
    //         length l of DAC, Kummer surface K 
    // Outputs: Dset[2][0] + [m]*Dset[1][0] + [n]*Dset[1][1]
    ////////////////////////////////////////////////////////////////////

    
    P,D,DT := Explode(Dset);
    first_add,b0,b1,b2,b3:=ENCODE(a,b,len); // compute 3-dim addition chain a la DJB
                        
    I:=[[1,1], [2,2], [1,1]];                   

    if first_add eq 1 then          
        P[3]:=PseudoAddKummer(P[3],D[2],D[1],K); I[3][2]+:=1;
    else
        P[3]:=PseudoAddKummer(P[3],D[1],D[2],K); I[3][1]+:=1;
    end if;

    for i:=len-1 to 1 by -1 do

        case [b0[i],b1[i],b2[i],b3[i]]:

            when [0,0,0,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[2][1],2*I[2][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[2],P[1],D[3],P[3],P[2],D[2],P[2],P[4],DT[IND(I)],K);
        
            when [0,0,0,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[2][1],2*I[2][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[2],P[1],D[3],P[3],P[2],D[1],P[2],P[4],DT[IND(I)],K);

            when [0,0,1,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[2][1],2*I[2][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[2],P[1],D[4],P[3],P[2],D[2],P[2],P[4],DT[IND(I)],K);

            when [0,0,1,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[2][1],2*I[2][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[2],P[1],D[4],P[3],P[2],D[1],P[2],P[4],DT[IND(I)],K);

            when [0,1,0,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[1][1],2*I[1][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[1],P[2],D[3],P[3],P[1],D[2],P[1],P[4],DT[IND(I)],K);

            when [0,1,0,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[1][1],2*I[1][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[1],P[2],D[3],P[3],P[1],D[1],P[1],P[4],DT[IND(I)],K);

            when [0,1,1,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[1][1],2*I[1][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[1],P[2],D[4],P[3],P[1],D[2],P[1],P[4],DT[IND(I)],K);

            when [0,1,1,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[1][1],2*I[1][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[1],P[3],P[4]:=DoubleThriceAdd(P[1],P[2],D[4],P[3],P[1],D[1],P[1],P[4],DT[IND(I)],K);

            when [1,0,0,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[2],D[2],P[1],P[2],D[3],P[3],P[4],DT[IND(I)],K);

            when [1,0,0,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[1],D[1],P[1],P[2],D[3],P[3],P[4],DT[IND(I)],K);

            when [1,0,1,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[2],D[2],P[1],P[2],D[4],P[3],P[4],DT[IND(I)],K);

            when [1,0,1,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[1],D[1],P[1],P[2],D[4],P[3],P[4],DT[IND(I)],K);

            when [1,1,0,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[1],D[2],P[1],P[2],D[3],P[3],P[4],DT[IND(I)],K);

            when [1,1,0,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[2],D[1],P[1],P[2],D[3],P[3],P[4],DT[IND(I)],K);

            when [1,1,1,0]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[1][1]+I[3][1],I[1][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[1],D[2],P[1],P[2],D[4],P[3],P[4],DT[IND(I)],K);

            when [1,1,1,1]: I:=[[I[1][1]+I[2][1],I[1][2]+I[2][2]], [2*I[3][1],2*I[3][2]], [I[2][1]+I[3][1],I[2][2]+I[3][2]]];
            P[2],P[3],P[1],P[4]:=DoubleThriceAdd(P[3],P[2],D[1],P[1],P[2],D[4],P[3],P[4],DT[IND(I)],K);

        end case;

    end for;

    return P[4];    

end function;
