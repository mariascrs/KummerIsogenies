clear;

// Here you can input your parameters 
e2:=60;
e3:=40;
f:=7;
p:=f*2^e2*3^e3-1;
Fp := GF(p);
Fp2<i>:=ExtensionField<Fp,x|x^2+1>;

poly<X,Y,Z,T> := PolynomialRing(Fp2,4);

/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

load "functions.m";

// Here we set up the Kummer surfaces and the 3-torsion points R, S 

_<x> := PolynomialRing(Fp2);
C:=HyperellipticCurve(x^6+1);
C:=random_walk(C);
K,rosen:=KummerFromCurve(C);
C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
J:=Jacobian(C);

// Getting the points of order 3 (we do this on the Jacobian and then push it down to the Kummer surface)
// We want these points to generate the kernel of a (3,3)-isogeny so we also require that their Weil pairing is trivial
// (i.e. generate a maximal isotropic group)
repeat
    repeat
        R:=((p+1) div 3)*Random(J);
    until 3*R eq J!0 and R ne J!0;
    repeat
        S:=((p+1) div 3)*Random(J);
    until 3*S eq J!0 and S notin [J!0,R,-R];
until WeilPairing(R,S,3) eq 1;

// Push the points R, S to the Kummer Surface K (rosen are the Rosenhein Invariants of the corresponding curve C)
R:=JtoK(R,rosen,K);
S:=JtoK(S,rosen,K);


/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

load "iso33.m";
"";

thetas:=K[2];

P:=[X,Y,Z,T];

a,b,c,d := Explode(thetas);
XR,YR,ZR,TR := Explode(R);
XS,YS,ZS,TS := Explode(S);


P3 := ProjectiveSpace(Fp2, 3);
Kum := Surface(P3, ComputeKummer(K));

"Getting formulae for intersection where the coefficients are symmetric in R, S...\n";
"psi: ";
psi := iso33sym(R,S,P,K);
normalise_phi(psi);
"Are the coordinates of psi symmetric?";
is_iso33sym_symmetric(R,S,P,K);
"";

psi := map<Kum -> P3 | psi>;
Kum1<X,Y,Z,T> := Image(psi);
Kum1:=DefiningPolynomials(Kum1)[1];
image_thetas := psi(thetas);
image_thetas := [image_thetas[i] : i in [1..4]];
"The image is: ";
Kum1;
"Is this in fast Kummer surface form?";
check := CheckFastKummerForm(poly!Kum1, image_thetas);
check;
"";

if not check then
    "Let's put this in fast Kummer surface form by the linear transformation (X,Y,Z,T) --> ( X/alpha1, 2*Y/alpha2, 2*Z/alpha1, 2*T/(3*alpha2) )...\n";
    "phi: ";

    phi33 := iso33(R,S,P,K);
    normalise_phi(phi33);
    
    phi33 := map<Kum -> P3 | phi33>;
    Kum2<X,Y,Z,T> := Image(phi33);
    Kum2:=DefiningPolynomials(Kum2)[1];
    image_thetas := phi33(thetas);
    image_thetas := [image_thetas[i] : i in [1..4]];
    "The image is: ";
    Kum2;
    "Is this in Fast Kummer surface form?";
    CheckFastKummerForm(poly!Kum2, image_thetas);
    "";
end if;
