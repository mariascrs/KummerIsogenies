clear;

/////////////////////////////////////
//////// Loading functions //////////
/////////////////////////////////////

load "kummer_arithmetic.m";
load "isogeny33.m";


/////////////////////////////////////////////
////// Hash without optimal strategies //////
/////////////////////////////////////////////

hash:=function(data,gens,TC,K)

    a,b,c,len,e:=Explode(data);
    P1,D1,DT1,P2,D2,DT2:=Explode(gens);

    //Q1:=P1+a*P3+b*P4; 
    Q1:=threeDAC([P1,D1,DT1],a,b,len,K);
    //Q2:=P2+b*P3+c*P4;
    Q2:=threeDAC([P2,D2,DT2],b,c,len,K);

    for ee:=e-1 to 1 by -1 do

        R:=Q1; S:=Q2;

        // "e: " cat IntegerToString(ee);

        for tri:=1 to ee do
            R:=TripleKummer(R,TC);
            S:=TripleKummer(S,TC);
        end for;

        phiPts,image_thetas:=Isogeny33ImageAndPoints(R,S,TC,[[Q1,Q2]]);

        Q1,Q2:=Explode(&cat phiPts);
        TC:=TriplingConstantsFromThetas(image_thetas);

    end for;

    return image_thetas;

end function;


/////////////////////////////////////////////
/////// Hash with optimal strategies ////////
/////////////////////////////////////////////

hash_optimal:=function(data,gens,TC,K,strategy)

    a,b,c,len,e:=Explode(data);
    P1,D1,DT1,P2,D2,DT2:=Explode(gens);

    //Q1:=P1+a*P3+b*P4; 
    Q1:=threeDAC([P1,D1,DT1],a,b,len,K);
    //Q2:=P2+b*P3+c*P4;
    Q2:=threeDAC([P2,D2,DT2],b,c,len,K);

    pts:=[];
    index:=0;
    indices:=[];
    R:=Q1; S:=Q2;

    for row:=1 to e-1 do

        while index lt (e-row) do
            Append(~pts,[R,S]);
            Append(~indices,index);
            m := strategy[e-index-row+1];

            for tri:=1 to m do
                R:=TripleKummer(R,TC);
                S:=TripleKummer(S,TC);
            end for;

            index +:= m;
        end while;

        pts,image_thetas:=Isogeny33ImageAndPoints(R,S,TC,pts);
        TC:=TriplingConstantsFromThetas(image_thetas);

        if #pts gt 0 then
            assert #indices eq #pts;
            R:=pts[#pts][1];
            S:=pts[#pts][2];
            index:=indices[#pts];        
            Prune(~pts);
            Prune(~indices);
        end if;

    end for;

    return image_thetas;

end function;

/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////


load "params128.m";
//load "params192.m";
//load "params256.m";

l:=Ceiling(Log(2,3^e));
len:=l+1;

times := [];
for i in [1..100] do
	"Sample: ", i;
	a:=2^l+Random(1,3^(e-1)); 
    b:=2^l+Random(1,3^(e-1));
    c:=2^l+Random(1,3^(e-1));
    data:=[a,b,c,len,e];
    TC:=TriplingConstantsFromThetas(K[2]);
	t:=Cputime();
    // h1 := hash(data,gens,TC,K);
	h2 := hash_optimal(data,gens,TC,K,strategy);
    // normalise(h1) eq normalise(h2);
	Append(~times, Cputime(t));
end for;

"Average:", &+times/100;
