Isogeny33Evaluate:=function(P,u)

    /////////////////////////////////////////////////////////////////
    // Computes the image of a point under a (3,3)-isogeny with coefficients [u1,u2,u3,u4,u5]
    // Inputs: a tuple of four projective field elements [x,y,z,t]
    // Outputs: Coefficients coeffs = [u1,u2,u3,u4,u5]
    /////////////////////////////////////////////////////////////////


    u1,u2,u3,u4,u5:=Explode(u);
    X,Y,Z,T:=Explode(P);
    X2,Y2,Z2,T2:=Explode(Squaring([X,Y,Z,T]));                             

    u5XY:=u5*X*Y;                                                           
    u5ZT:=u5*Z*T;                                                           

    p1:= X*(u1*X2+u2*Y2+u3*Z2+u4*T2)+u5ZT*Y;                                
    p2:= Y*(u1*Y2+u2*X2+u3*T2+u4*Z2)+u5ZT*X;                                
    p3:= Z*(u1*Z2+u2*T2+u3*X2+u4*Y2)+u5XY*T;                                
    p4:= T*(u1*T2+u2*Z2+u3*Y2+u4*X2)+u5XY*Z;                                

    return [p1,p2,p3,p4];

end function;

Isogeny33ImageAndPoints := function(R,S,TC,pts)

    ////////////////////////////////////////////////////////////////////
    // Given generators of a (3,3)-subgroup computes the 
    // corresponding isogeny, and pushes any points through the isogeny.
    // Inputs: Generators <R,S> of (3,3)-subgroup, Tripling constants TC 
    //         needed to comput isogeny, Points pts, given as a tuple of 
    //         four projective field elements [x,y,z,t], to be pushed 
    //         through the isogeny
    // Outputs: The image of the points in pts under the (3,3)-isogeny,
    //          the thetas constants of the image Kummer surface 
    ////////////////////////////////////////////////////////////////////
    
    a:=TC[1][1]; b:=TC[1][2]; c:=TC[1][3]; d:=TC[1][4];
    XR:=R[1]; YR:=R[2]; ZR:=R[3]; TR:=R[4];
    XS:=S[1]; YS:=S[2]; ZS:=S[3]; TS:=S[4];

    sO:=TC[2];
    hO:=TC[3];                              
    hOi:=TC[5];                                      

    ab:=a*b; ac:=a*c; ad:=a*d; bc:=b*c; bd:=b*d; cd:=c*d;           
    
    D12:=(ab-cd)*(ab+cd)*(ac-bd)*(ac+bd);                           
    D13:=(ab-cd)*(ab+cd)*(ad-bc)*(ad+bc);                           
    D23:=(ac-bd)*(ac+bd)*(ad-bc)*(ad+bc);                           

    hR:=Hadamard(Squaring(R));                                              
    hS:=Hadamard(Squaring(S));                                              
    hRO:=FourWayMult(hR,hOi);                                               
    hSO:=FourWayMult(hS,hOi);                                               
    HR:=Hadamard(hRO); 
    HS:=Hadamard(hSO); 
    RXX:=HR[1]; RYY:=HR[2]; RZZ:=HR[3]; RTT:=HR[4];
    SXX:=HS[1]; SYY:=HS[2]; SZZ:=HS[3]; STT:=HS[4];

    XYR:=XR*YR; XZR:=XR*ZR; XTR:=XR*TR; YZR:=YR*ZR; YTR:=YR*TR; ZTR:=ZR*TR; 
    XYS:=XS*YS; XZS:=XS*ZS; XTS:=XS*TS; YZS:=YS*ZS; YTS:=YS*TS; ZTS:=ZS*TS; 
       
    Rsig:=(ab-cd)*(ZTR+XYR);                                                            
    Rdel:=(ab+cd)*(XYR-ZTR);                                                            
    Ssig:=(ab-cd)*(ZTS+XYS);                                                            
    Sdel:=(ab+cd)*(XYS-ZTS);                                                            

    RXY:=(Rsig+Rdel); RZT:=(Rsig-Rdel); SXY:=(Ssig+Sdel); SZT:=(Ssig-Sdel); 

    Rsig:=(ac-bd)*(YTR+XZR);                                                            
    Rdel:=(ac+bd)*(XZR-YTR);                                                            
    Ssig:=(ac-bd)*(YTS+XZS);                                                            
    Sdel:=(ac+bd)*(XZS-YTS);                                                            

    RXZ:=(Rsig+Rdel);  RYT:=(Rsig-Rdel); SXZ:=(Ssig+Sdel); SYT:=(Ssig-Sdel);

    Rsig:=(ad-bc)*(XTR+YZR);                                                            
    Rdel:=(ad+bc)*(XTR-YZR);                                                            
    Ssig:=(ad-bc)*(XTS+YZS);                                                            
    Sdel:=(ad+bc)*(XTS-YZS);                                                            

    RXT:=(Rsig+Rdel); RYZ:=(Rsig-Rdel); SXT:=(Ssig+Sdel); SYZ:=(Ssig-Sdel); 

    c2:=D23*SXY+D13*SXZ+D12*SXT;                                                
    d2:=D23*RXY+D13*RXZ+D12*RXT;                                                

    E1:=(d2*SZT-c2*RZT)*D23;                                                
    E2:=SXX*RYY-RXX*SYY;                                                   

    E1R:=E1*RXX;                                                            
    E1S:=E1*SXX;                                                           
    E2c:=E2*c2;                                                             
    E2d:=E2*d2;                                                            

    u1:=E1R*SXX;                                                             
    u1:=u1+u1; 
    u2:=E1S*RYY+E1R*SYY+(E2d*SZT+E2c*RZT)*D23;                              
    u3:=E1S*RZZ+E1R*SZZ+(E2d*SYT+E2c*RYT)*D13;                              
    u4:=E1S*RTT+E1R*STT+(E2d*SYZ+E2c*RYZ)*D12;                              
    u5:=E2c*d2;                                                             
    u5:=u5+u5;

    a2,b2,c2,d2:=Explode(sO);                              

    u5XY:=u5*ab;                                                           
    u5ZT:=u5*cd;                                                           

    aa:= a*(u1*a2+u2*b2+u3*c2+u4*d2)+u5ZT*b;                                
    bb:= b*(u1*b2+u2*a2+u3*d2+u4*c2)+u5ZT*a;                                
    cc:= c*(u1*c2+u2*d2+u3*a2+u4*b2)+u5XY*d;                                
    dd:= d*(u1*d2+u2*c2+u3*b2+u4*a2)+u5XY*c;                                

    u:=[u1,u2,u3,u4,u5];

    phiPts:=[[Isogeny33Evaluate(p,u) : p in P ] : P in pts];
    
    return phiPts,[aa,bb,cc,dd];

end function; 