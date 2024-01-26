//// Basic operations on the Kummer surface

Hadamard:=function(P);

    // Given a point P on Kummer surface K, outputs the Hadamard of this point 

    T1:=P[1]+P[2]; T2:=P[3]+P[4]; T3:=P[1]-P[2]; T4:=P[3]-P[4];

    return [T1+T2,T1-T2,T3+T4,T3-T4];

end function;

Squaring:=function(P);

    // Given a point P, outputs the coordinate-wise squaring 

    return [P[1]^2,P[2]^2,P[3]^2,P[4]^2];

end function;

FourWayMult:=function(V,W)

    // Given points V and W, outputs their coordinate-wise multiplication 

    return [V[1]*W[1],V[2]*W[2],V[3]*W[3],V[4]*W[4]];

end function;

Invert4Constants:=function(P) // 6M

    // Given a point P on Kummer surfac K, returns projective equivalent to [1/P[1],1/P[2],1/P[3],1/P[4]]
    

    pi1:=P[3]*P[4];
    pi2:=pi1*P[1];
    pi1:=pi1*P[2];
    pi3:=P[1]*P[2];
    pi4:=pi3*P[3];
    pi3:=pi3*P[4];

    return [pi1,pi2,pi3,pi4];

end function;



normalise:=function(P)

    // Given a point P = [P1 : P2 : P3 : P4], outputs the normalised point

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

ComputeKummer:=function(K);
    // Given Kummer surface K, computes the equation of the Kummer surface
	P3<X,Y,Z,T> := ProjectiveSpace(Fp2, 3);
	E:=K[1][1]; F:=K[1][2]; G:=K[1][3]; H:=K[1][4];

	return ((X^4+Y^4+Z^4+T^4)+ 2*E*X*Y*Z*T - F*(X^2*T^2+Y^2*Z^2)-G*(X^2*Z^2+Y^2*T^2) - H*(X^2*Y^2+T^2*Z^2));
end function;

//// Operations related to computation of (3,3)-isogenies

TriplingConstantsFromThetas:=function(O) 

    // Given theta constants, outputs the tripling constants needed for triping and computing (3,3)-isogeny formulae
    
    SO:=Squaring(O); 
    HSO:=Hadamard(SO); 
    K3:=Invert4Constants(O);
    K4:=Invert4Constants(HSO);

    return [O,SO,HSO,K3,K4];

end function;

normalise_phi := function(phi)
    // Given an isogeny phi = (phi1 : phi2 : phi3 : phi4), computes the normalised isogeny formulae (by dividing by first coefficient of phi1)
    p1,p2,p3,p4 := Explode(phi);

    c := Coefficients(p1)[1];
    return [p1/c, p2/c, p3/c, p4/c];
end function;

CheckFastKummerForm:=function(Kum_eqn, thetas)
	
    // Checks that a Kummer surface with equation Kum_eqn is in fast Kummer surface form (as given in Equation (1) in the accompanying paper)

    poly<X,Y,Z,T> := Parent(Kum_eqn);
	
    
	c1 := MonomialCoefficient(Kum_eqn,X*Y*Z*T);
	c2 := MonomialCoefficient(Kum_eqn,X^2*T^2);
	c3 := MonomialCoefficient(Kum_eqn,X^2*Z^2);
	c4 := MonomialCoefficient(Kum_eqn,X^2*Y^2);

	bool1 := c2 eq MonomialCoefficient(Kum_eqn,Y^2*Z^2);
	bool2 :=  c3 eq MonomialCoefficient(Kum_eqn,Y^2*T^2);
	bool3 := c4 eq MonomialCoefficient(Kum_eqn,Z^2*T^2);

	E:= c1/2;
	F:= -c2;
	G:= -c3;
	H:= -c4;

	a,b,c,d:=Explode(thetas);

	E1 := a*b*c*d*(a^2+b^2+c^2+d^2)*(d^2-a^2+b^2-c^2)*(d^2-a^2-b^2+c^2)*(d^2+a^2-b^2-c^2)/(a^2*d^2-b^2*c^2)/(b^2*d^2-a^2*c^2)/(c^2*d^2-a^2*b^2);
	F1:=(a^4+d^4 - b^4-c^4)/(a^2*d^2-b^2*c^2);
	G1:=(b^4+d^4 - a^4-c^4)/(b^2*d^2-a^2*c^2);
	H1:=(c^4+d^4 - a^4-b^4)/(c^2*d^2-a^2*b^2);

	
    bool4 := [E,F,G,H] eq [E1, F1, G1, H1];

	return bool1 and bool2 and bool3 and bool4;

end function;

//// Functions needed for setup

random_richelot := function(C)
    // Given a genus-2 curve, returns a pseudorandom genus-2 curve CC, where Jac(C) and Jac(CC) are (2,2)-isogenous
    isos:=RichelotIsogenousSurfaces(C);
	repeat
		CC:=isos[Random(1,#isos)];
	until Type(CC) eq CrvHyp;
	return CC;
end function;

random_walk := function(C)
    // Given a genus-2 curve C, ouputs a pseudorandom genus-2 curve CC such that Jac(C) and Jac(CC) are connected by a (2^r, 2^r)-isogeny, where r is a random integer between 5 and 12

	steps:=Random(5,12);
	CC:=C;
	for i:=1 to steps do
		CC:=random_richelot(CC);
	end for;
	return CC;
end function;



all_rosenhains:=function(C)
    // Given a genus-2 curve C, outputs all possible Rosenhain invariants

    f:=HyperellipticPolynomials(C);
    poly<x>:=Parent(f);
    assert Degree(f) eq 6;
    rts:=Roots(f); rts:=[rts[i][1]: i in [1..6]];
    set:={};

    for r1 in rts do 
        others1:=Remove(rts,Index(rts,r1));
        for r2 in others1 do
            others2:=Remove(others1,Index(others1,r2));
            for r3 in others2 do
                others3:=Remove(others2,Index(others2,r3));
                g:=(x-r1)/(x-r2)*(r3-r2)/(r3-r1);
                for r4 in others3 do
                    others4:=Remove(others3,Index(others3,r4));
                    lambda:=Evaluate(g,r4);
                    for r5 in others4 do
                        mu:=Evaluate(g,r5);
                        others5:=Remove(others4,Index(others4,r5));
                        nu:=Evaluate(g,others5[1]);
                        Include(~set,[lambda,mu,nu]);
                    end for;
                end for;
            end for;
        end for;
    end for;

    return set;

end function;

one_squared_thetas:=function(C)

    // Given a genus-2 curve C, outputs the Rosenhain invariants and the corresponding squared theta constants

    rosen:=all_rosenhains(C);

    poly<x>:=Parent(HyperellipticPolynomials(C));

    thetas:={};

    repeat 

		lmn:=Random(rosen);
		//lmn:=[Fp2!-1,Fp2!-3,Fp2!3];
        lambda:=lmn[1]; mu:=lmn[2]; nu:=lmn[3];
        CC:=HyperellipticCurve(x*(x-1)*(x-lambda)*(x-mu)*(x-nu));
        assert IsSquare(lambda*mu/nu);
        assert IsSquare(mu*(mu-1)*(lambda-nu)/(nu*(nu-1)*(lambda-mu)));

        d2:=1;
		c2sign:=Random([-1,1]);
		b2sign:=Random([-1,1]);
		c2:=c2sign*Sqrt(lambda*mu/nu);
        b2:=b2sign*Sqrt(mu*(mu-1)*(lambda-nu)/(nu*(nu-1)*(lambda-mu)));
        a2:=b2*c2*nu/mu;

    until IsSquare(a2) and IsSquare(b2) and IsSquare(c2);

    return [a2,b2,c2,d2],[lambda,mu,nu];

end function;

KummerFromCurve:=function(C)
    // Given a genus-2 curve, output the corresponding fast Kummer surface

	sq_thetas,rosen:=one_squared_thetas(C);
	
	a:=Sqrt(sq_thetas[1]);
	b:=Sqrt(sq_thetas[2]);
	c:=Sqrt(sq_thetas[3]);
	d:=Sqrt(sq_thetas[4]);

	E := a*b*c*d*(a^2+b^2+c^2+d^2)*(d^2-a^2+b^2-c^2)*(d^2-a^2-b^2+c^2)*(d^2+a^2-b^2-c^2)/(a^2*d^2-b^2*c^2)/(b^2*d^2-a^2*c^2)/(c^2*d^2-a^2*b^2);
	F:=(a^4+d^4 - b^4-c^4)/(a^2*d^2-b^2*c^2);
	G:=(b^4+d^4 - a^4-c^4)/(b^2*d^2-a^2*c^2);
	H:=(c^4+d^4 - a^4-b^4)/(c^2*d^2-a^2*b^2);

	K1:=[E,F,G,H];
	K2:=[a,b,c,d];
	K3:=Invert4Constants(K2);
	K4:=Invert4Constants(Hadamard(sq_thetas));

	K:=[K1,K2,K3,K4];

	return K,rosen;

end function;


JtoK:=function(P,rosen,K)

	// Given a point on the Jacobian J, outputs (2 times) the corresponding point on the Kummer surface

	a2:=K[2][1]^2; b2:=K[2][2]^2; c2:=K[2][3]^2; d2:=K[2][4]^2; 
	lambda:=rosen[1]; mu:=rosen[2]; nu:=rosen[3];

	q:=Coefficients(P[1])[2];
	r:=Coefficients(P[1])[1];
	t:=Coefficients(P[2])[1];

	X:=a2*(r*(mu-r)*(lambda+q+nu)-t^2);
	Y:=b2*(r*(nu*lambda-r)*(1+q+mu)-t^2);
	Z:=c2*(r*(nu-r)*(lambda+q+mu)-t^2);
	T:=d2*(r*(mu*lambda-r)*(1+q+nu)-t^2);

	P:=[X,Y,Z,T];

	P:=Hadamard(P);
	P:=Squaring(P);
	P:=FourWayMult(P,K[4]);
	P:=Hadamard(P);
	P:=FourWayMult(P,K[3]);

	return normalise(P);

end function;