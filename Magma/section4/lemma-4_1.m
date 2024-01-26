// This file aims to verify the clain in Lemma 4.1 that, from the equations arising from Equation (5), we obtain the form of the isogeny in the statement of the lemma


// Setting up the our polynomial ring containing all the coefficients of the cubic monomials in phi
pp<a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,
    b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,
    c01,c02,c03,c04,c05,c06,c07,c08,c09,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
    d01,d02,d03,d04,d05,d06,d07,d08,d09,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20> := PolynomialRing(Rationals(), 20*4);

coeffs := [
            [a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20], 
            [b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20],
            [c01,c02,c03,c04,c05,c06,c07,c08,c09,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20],
            [d01,d02,d03,d04,d05,d06,d07,d08,d09,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20]
        ];

// Setting polynomial ring for isogeny formulae
poly<X,Y,Z,T> := PolynomialRing(pp, 4);

// Setting up \sigma_i(X:Y:Z:T) where \sigma_i is the action of the 2 torsion point T_i (as defined in Section 2 of the paper)
coord_translates := [
    [X,Y,Z,T],[X,Y,-Z,-T],[X,-Y,Z,-T],[X,-Y,-Z,T],
    [Y,X,T,Z],[Y,X,-T,-Z], [Y,-X,T,-Z],[Y,-X,-T,Z],
    [Z,T,X,Y],[Z,T,-X,-Y],[Z,-T,X,-Y],[Z,-T,-X,Y],
    [T,Z,Y,X],[T,Z,-Y,-X],[T,-Z,Y,-X],[T,-Z,-Y,X]];


// Setting up the monomials of degree 3
mons := MonomialsOfDegree(poly, 3);

"Setting up phi = [phi1, phi2, phi3, phi4]:\n";
phi1 := poly!&+[coeffs[1][i]*mons[i] : i in [1..20]];
phi2 := poly!&+[coeffs[2][i]*mons[i] : i in [1..20]];
phi3 := poly!&+[coeffs[3][i]*mons[i] : i in [1..20]];
phi4 := poly!&+[coeffs[4][i]*mons[i] : i in [1..20]];

// Setting up \sigma_i((phi1 : phi2 : phi3 : phi4))

phi_translates := [
    [phi1,phi2,phi3,phi4],[phi1,phi2,-phi3,-phi4],[phi1,-phi2,phi3,-phi4],[phi1,-phi2,-phi3,phi4],
    [phi2,phi1,phi4,phi3],[phi2,phi1,-phi4,-phi3],[phi2,-phi1,phi4,-phi3],[phi2,-phi1,-phi4,phi3],
    [phi3,phi4,phi1,phi2],[phi3,phi4,-phi1,-phi2],[phi3,-phi4,phi1,-phi2],[phi3,-phi4,-phi1,phi2],
    [phi4,phi3,phi2,phi1],[phi4,phi3,-phi2,-phi1],[phi4,-phi3,phi2,-phi1],[phi4,-phi3,-phi2,phi1]
    ];

phi := phi_translates[1];

"phi1: ", phi1;
"phi2: ", phi2;
"phi3: ", phi3;
"phi4: ", phi4;
"";


// Finding relations between the coefficients in the isogeny formulae

"Finding relations..\n";
rels := [];

for p in [1..4] do
    for i in [2..16] do
        F := Evaluate(phi[p], coord_translates[i]) - phi_translates[i][p];
        R := [MonomialCoefficient(F, m) : m in mons];
        for r in R do
            if r ne 0 and r notin rels then 
                Append(~rels, r);
            end if;
        end for;
    end for;
end for;


// Using these relations to reduce the number of coefficients
cs := &cat coeffs;
M:=[];
for r in rels do
    Append(~M, [MonomialCoefficient(r, c) : c in cs]);
end for;
M:=Matrix(M);
kermat := Basis(Nullspace(Transpose(M)));

new_cs := [pp!0 : i in [1..80]];
for k in kermat do 
    i:=0;
    tmp := 0;
    while tmp eq 0 do
        i +:= 1;
        tmp := k[i];
    end while;
    for j in [1..80] do
        if k[j] ne 0 then 
            new_cs[j] := cs[i];
        end if;
    end for;
end for;

// Outputting the new coordinates of the isogeny

phi1 := poly!&+[new_cs[i]*mons[i] : i in [1..20]];
phi2 := poly!&+[new_cs[20+i]*mons[i] : i in [1..20]];
phi3 := poly!&+[new_cs[40+i]*mons[i] : i in [1..20]];
phi4 := poly!&+[new_cs[60+i]*mons[i] : i in [1..20]];

"The coordinates of phi are now";


"phi1: ", phi1;
"phi2: ", phi2;
"phi3: ", phi3;
"phi4: ", phi4;
"";
