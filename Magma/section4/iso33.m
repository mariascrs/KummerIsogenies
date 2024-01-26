
biquadratics:=function(thetas,P,Q)

    // Given theta constants of a Kummer surface K, and points P, Q on K, outputs the biquadratics B_{ij}(P,Q)

    a,b,c,d := Explode(thetas);
    a:=K[2][1]; b:=K[2][2]; c:=K[2][3]; d:=K[2][4];
    XP:=P[1]; YP:=P[2]; ZP:=P[3]; TP:=P[4];
    XQ:=Q[1]; YQ:=Q[2]; ZQ:=Q[3]; TQ:=Q[4];  

    B11:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))+(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))+(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))+(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B22:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))+(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))-(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))-(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B33:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))-(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))+(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))-(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B44:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))-(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))-(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))+(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));

    B12:=(2*(a*b*(TP*TQ*ZP*ZQ+XP*XQ*YP*YQ)-c*d*(TP*XQ*YQ*ZP+TQ*XP*YP*ZQ)))/(((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)-((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B13:=(2*(a*c*(TP*TQ*YP*YQ+XP*XQ*ZP*ZQ)-b*d*(TP*XQ*YP*ZQ+TQ*XP*YQ*ZP)))/(((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B14:=(2*(a*d*(TP*TQ*XP*XQ+YP*YQ*ZP*ZQ)-b*c*(TP*XP*YQ*ZQ+TQ*XQ*YP*ZP)))/(((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2));

    B23:=(2*(b*c*(TP*TQ*XP*XQ+YP*YQ*ZP*ZQ)-a*d*(TP*XP*YQ*ZQ+TQ*XQ*YP*ZP)))/(((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B24:=(2*(b*d*(TP*TQ*YP*YQ+XP*XQ*ZP*ZQ)-a*c*(TP*XQ*YP*ZQ+TQ*XP*YQ*ZP)))/(((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2));

    B34:=(2*(c*d*(TP*TQ*ZP*ZQ+XP*XQ*YP*YQ)-a*b*(TP*XQ*YQ*ZP+TQ*XP*YP*ZQ)))/(((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2));

    B21:=B12;
    B31:=B13;
    B41:=B14;

    B32:=B23;
    B42:=B24;

    B43:=B34;

    B:=[
        [B11,B12,B13,B14],
        [B21,B22,B23,B24],
        [B31,B32,B33,B34],
        [B41,B42,B43,B44]
       ];

    return B;

end function;


iso33sym:=function(R,S,P, K)

    // Given a Kummer surface K, 3-torsion points R,S on K generating a (3,3)-subgroup, and a generic point P = [X:Y:Z:T], outputs formulae for the intersection (as described in Section 4 of the accompanying paper) where the coordinates are symmetric.

    thetas := K[2];
    a,b,c,d := Explode(thetas);
    XR,YR,ZR,TR := Explode(R);
    XS,YS,ZS,TS := Explode(S);

    R_mat:=biquadratics(thetas,P,R);
    S_mat:=biquadratics(thetas,P,S);

    hR:=Hadamard(Squaring(R));
    hS:=Hadamard(Squaring(S));
    hO:=Hadamard(Squaring(K[2]));
    hRO:=[hR[i]/hO[i]: i in [1..4]];
    hSO:=[hS[i]/hO[i]: i in [1..4]];
    HR:=Hadamard(hRO);
    HS:=Hadamard(hSO);

    RXX:=Coefficient(R_mat[1,1],X,2);
    RYY:=Coefficient(R_mat[1,1],Y,2);
    RZZ:=Coefficient(R_mat[1,1],Z,2);
    RTT:=Coefficient(R_mat[1,1],T,2);
    RXY:=Coefficient(Coefficient(R_mat[1,2],X,1),Y,1);
    RZT:=Coefficient(Coefficient(R_mat[1,2],Z,1),T,1);
    RXZ:=Coefficient(Coefficient(R_mat[1,3],X,1),Z,1);
    RYT:=Coefficient(Coefficient(R_mat[1,3],Y,1),T,1);
    RXT:=Coefficient(Coefficient(R_mat[1,4],X,1),T,1);
    RYZ:=Coefficient(Coefficient(R_mat[1,4],Y,1),Z,1);

    SXX:=Coefficient(S_mat[1,1],X,2);
    SYY:=Coefficient(S_mat[1,1],Y,2);
    SZZ:=Coefficient(S_mat[1,1],Z,2);
    STT:=Coefficient(S_mat[1,1],T,2);
    SXY:=Coefficient(Coefficient(S_mat[1,2],X,1),Y,1);
    SZT:=Coefficient(Coefficient(S_mat[1,2],Z,1),T,1);
    SXZ:=Coefficient(Coefficient(S_mat[1,3],X,1),Z,1);
    SYT:=Coefficient(Coefficient(S_mat[1,3],Y,1),T,1);
    SXT:=Coefficient(Coefficient(S_mat[1,4],X,1),T,1);
    SYZ:=Coefficient(Coefficient(S_mat[1,4],Y,1),Z,1);


    c2:=SXY+SXZ+SXT;
    d2:=RXY+RXZ+RXT;
   
    p1XXX:=(d2*SZT-c2*RZT)*(HS[1]*RXX+HR[1]*SXX);
    p1XYY:=(d2*SZT-c2*RZT)*(HS[1]*RYY+HR[1]*SYY)+(HS[1]*RYY-HR[1]*SYY)*(c2*RZT+d2*SZT);
    p1XZZ:=(d2*SZT-c2*RZT)*(HS[1]*RZZ+HR[1]*SZZ)+(HS[1]*RYY-HR[1]*SYY)*(c2*RYT+d2*SYT);
    p1XTT:=(d2*SZT-c2*RZT)*(HS[1]*RTT+HR[1]*STT)+(HS[1]*RYY-HR[1]*SYY)*(c2*RYZ+d2*SYZ);
    p1YZT:=                                     +(HS[1]*RYY-HR[1]*SYY)*2*c2*d2;
    
    phi1 := p1XXX*X^3+p1XYY*X*Y^2+p1XZZ*X*Z^2+p1XTT*X*T^2+p1YZT*Y*Z*T;

    c2:=SXY+SXZ+SXT;
    d2:=RXY+RXZ+RXT;

    p2YYY:=(d2*SYZ-c2*RYZ)*(SXX*RXX+RXX*SXX);
    p2YXX:=(d2*SYZ-c2*RYZ)*(SXX*RYY+RXX*SYY)+(SXX*RTT-STT*RXX)*(c2*RZT+d2*SZT);
    p2YZZ:=(d2*SYZ-c2*RYZ)*(SXX*RTT+RXX*STT)+(SXX*RTT-STT*RXX)*(c2*RYZ+d2*SYZ);
    p2YTT:=(d2*SYZ-c2*RYZ)*(SXX*RZZ+RXX*SZZ)+(SXX*RTT-STT*RXX)*(c2*RYT+d2*SYT);
    p2XZT:=                                 +(SXX*RTT-STT*RXX)*2*c2*d2;

    phi2 := p2YYY*Y^3+p2YXX*Y*X^2+p2YZZ*Y*Z^2+p2YTT*Y*T^2+p2XZT*X*Z*T;

    c2:=SXY+SXT+SXZ;
    d2:=RXY+RXT+RXZ;

    p3ZZZ:=(d2*SZT-c2*RZT)*(SXX*RXX+RXX*SXX);
    p3ZXX:=(d2*SZT-c2*RZT)*(SXX*RZZ+RXX*SZZ)+(SXX*RYY-RXX*SYY)*(c2*RYT+d2*SYT);
    p3ZYY:=(d2*SZT-c2*RZT)*(SXX*RTT+RXX*STT)+(SXX*RYY-RXX*SYY)*(c2*RYZ+d2*SYZ);
    p3ZTT:=(d2*SZT-c2*RZT)*(SXX*RYY+RXX*SYY)+(SXX*RYY-RXX*SYY)*(c2*RZT+d2*SZT);
    p3XYT:=                                 +(SXX*RYY-RXX*SYY)*2*c2*d2;

    phi3 := p3ZZZ*Z^3+p3ZXX*Z*X^2+p3ZYY*Z*Y^2+p3ZTT*Z*T^2+p3XYT*X*Y*T;

    c2:=SXY+SXT+SXZ;
    d2:=RXY+RXT+RXZ;

    p4TTT:=3*(d2*SYZ-c2*RYZ)*(SXX*RXX+RXX*SXX);
    p4TXX:=3*(d2*SYZ-c2*RYZ)*(SXX*RTT+RXX*STT)+3*(SXX*RTT-STT*RXX)*(c2*RYZ+d2*SYZ);
    p4TYY:=3*(d2*SYZ-c2*RYZ)*(SXX*RZZ+RXX*SZZ)+3*(SXX*RTT-STT*RXX)*(c2*RYT+d2*SYT);
    p4TZZ:=3*(d2*SYZ-c2*RYZ)*(SXX*RYY+RXX*SYY)+3*(SXX*RTT-STT*RXX)*(c2*RZT+d2*SZT);
    p4XYZ:=                                   +3*(SXX*RTT-STT*RXX)*2*c2*d2;

    phi4:= p4TTT*T^3+p4TXX*T*X^2+p4TYY*T*Y^2+p4TZZ*T*Z^2+p4XYZ*X*Y*Z;

    return [phi1,phi2,phi3,phi4];

end function;



iso33 := function(R,S,P, K)

    // Given a Kummer surface K, 3-torsion points R,S on K generating a (3,3)-subgroup, and a generic point P = [X:Y:Z:T], outputs formulae for the (3,3)-isogeny after scaling the intersection by M (as described in Section 4 of the paper)

    thetas := K[2];
    a,b,c,d := Explode(thetas);
    XR,YR,ZR,TR := Explode(R);
    XS,YS,ZS,TS := Explode(S);

    hO:=Hadamard(Squaring(thetas));
    D1:=hO[1]*hO[2]-hO[3]*hO[4];
    D2:=hO[1]*hO[3]-hO[2]*hO[4];
    D3:=hO[1]*hO[4]-hO[2]*hO[3];
    RXY:=XR*YR*a*b-TR*ZR*c*d; SXY:=XS*YS*a*b-TS*ZS*c*d;
    RZT:=TR*ZR*a*b-XR*YR*c*d; SZT:=TS*ZS*a*b-XS*YS*c*d;
    RXZ:=XR*ZR*a*c-TR*YR*b*d; SXZ:=XS*ZS*a*c-TS*YS*b*d;
    RYT:=TR*YR*a*c-XR*ZR*b*d; SYT:=TS*YS*a*c-XS*ZS*b*d;
    RXT:=TR*XR*a*d-YR*ZR*b*c; SXT:=TS*XS*a*d-YS*ZS*b*c;
    RYZ:=YR*ZR*a*d-TR*XR*b*c; SYZ:=YS*ZS*a*d-TS*XS*b*c;

    c2:=D2*D3*SXY+D1*D3*SXZ+D1*D2*SXT; 
    d2:=D2*D3*RXY+D1*D3*RXZ+D1*D2*RXT;


    alpha1 := (d2*SZT-c2*RZT)*D3;
    alpha2 := (d2*SYZ-c2*RYZ)*D1;
    psi := iso33sym(R,S,P,K);

    return [psi[1]/alpha1, 2*psi[2]/alpha2, 2*psi[3]/alpha1, 2*psi[4]/(3*alpha2)];

end function;               



is_iso33sym_symmetric:=function(R,S,P,K)

    // Checking that the intersection formulae given by (3,3)-subgroup <R, S> on Kummer surface K is symmetric

    return normalise_phi(iso33sym(R,S,P,K)) eq normalise_phi(iso33sym(S,R,P,K));

end function;


