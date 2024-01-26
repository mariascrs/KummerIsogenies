//clear;

//////////////////////////////////
// Parameters for secp = 192
//////////////////////////////////


e:=115; 
p:=37*2^4*3^e-1;
Fp:=GF(p);
Fp2<i>:=ExtensionField<Fp,x|x^2+1>;
_<x>:=PolynomialRing(Fp2);

//Top Strategy: [ 63, 52 ]
//375 TripleKummer and 417 Isogeny33Eval ==  27594 total units
strategy:=[ 0, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 11, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 18, 19, 20, 21, 21, 21, 21, 21, 22, 23, 24, 25, 26, 27, 27, 28, 29, 30, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 35, 36, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 48, 48, 48, 48, 48, 49, 50, 51, 52 ];

//load "setup.m";
//K,consts:=setup(f,4,e);

K:=[
    [
        2143365707512131625548481541494791671277397003451\
            568077225*i + 1640762889189872748133586452500\
            823238005062987004295576872,
        3384424124152531124702711559695453209869799956729\
            101623402*i + 1437045354495836548551312407804\
            458333531429472342861841685,
        5059910409857897425199825641246424352961128292297\
            84931194*i + 21975367870929247594650171130212\
            38056669162431481777651068,
        2644129362648307915216884727938480393346071397445\
            642360049*i + 1758246594949217824815417158277\
            750516740187653402822248764
    ],
    [
        2034077005136398522706550060631095799658714536589\
            89078292*i + 13940452183559784047836152995294\
            43979378271226503064092002,
        2426179325581650152269164583263966444791957641865\
            618676836*i + 2122894590735883871012489637721\
            179229225248713613291447957,
        3319983387200550564802614719676528310861486369861\
            451793591*i + 2709782916890541485848053028401\
            296204662613668923369965590,
        1
    ],
    [
        2325889772657614171076125196906630529710114340299\
            398297825*i + 5222987433491901126387838679819\
            11926380034545333897454108,
        3096448160155933464718565538688679472341281546692\
            380786340*i + 3934391540159745702180354022449\
            247809129290654095362346808,
        4090362466060298329788064980250940942659642699797\
            346575933*i + 1364778294567514477697798663557\
            085412446197170400668731471,
        2583397667594076309171840612013643184160203112910\
            554349837*i + 4047378812265433796080335688319\
            571815473655713231309458967
    ],
    [
        3845876715522380055813648378568559374301917345399\
            476905445*i + 1371994450579567697136495196094\
            437902986017715456989609337,
        1359841965389615333749449756421483530376259490040\
            899253046*i + 2795398025696002161088747253659\
            737449659800353207222759269,
        4965935450005665078971448854578587532006976310428\
            81884687*i + 19441644210491293169542146365696\
            20433313888166213356019912,
        7251603532257054695812924395728688184844630576519\
            08589434*i + 32641859382293888062504526655812\
            62647227880942626792901738
    ]
];

gens:=[
    [
        [
            305911177965252789833383994154493036996770951\
                0172993886444*i + 92565742612037761049150\
                4677137923505124932074435358013568,
            298080905845429690262707890078011680618909522\
                8301685298308*i + 34322372318479954444759\
                97009861549810901885286817560821549,
            995462354109678197105807353386077887495415751\
                245688653545*i + 192754851293894114748183\
                7522784332496852136496595867631368,
            1
        ],
        [
            322976807103332165146148181113248196660064932\
                5392559718602*i + 32128018232436803164244\
                60650119693628695876941007409008035,
            129506428974323366430162360304434146718789665\
                1732402736102*i + 93769068037729199848950\
                1472940615413224522039795090905695,
            228497068144326080559926525740325440139944192\
                2369921834090*i + 23981631512212292618880\
                92150759916000269326600272170787172,
            1
        ],
        [
            305911177965252789833383994154493036996770951\
                0172993886444*i + 92565742612037761049150\
                4677137923505124932074435358013568,
            298080905845429690262707890078011680618909522\
                8301685298308*i + 34322372318479954444759\
                97009861549810901885286817560821549,
            995462354109678197105807353386077887495415751\
                245688653545*i + 192754851293894114748183\
                7522784332496852136496595867631368,
            1
        ],
        [
            392654243114819572007316904624658772446650243\
                7600951398022*i + 14847595862242537578476\
                5106328673500086079810038810121091,
            191924663660827205986984400291025121557540597\
                8637141249506*i + 32276104511718154519641\
                26144152372126572647148655295434409,
            249954308951368844904502344583592741528408700\
                457262367342*i + 170393769018270342313115\
                9195400309544566611120177161803313,
            1
        ]
    ],
    [
        [
            318061222727715218215094343579627702176550032\
                89375539139*i + 1020719658316342498347286\
                503110742505967417215854384643316,
            313011583926264772514846877546187774690030488\
                3925282354478*i + 50511150760349870752306\
                3753765033348657300943396885453085,
            422057279207623869511260778469104421613550943\
                2729031297689*i + 21052163327757514892731\
                5881629447112886313269875222127189,
            1
        ],
        [
            254189867054532339101846669674285727138806408\
                6970879515442*i + 10583126258379777096348\
                38840884875080817693628368747427569,
            218759986458508666349438519622066986324977456\
                810349181996*i + 333545735753950379126666\
                2593941254062416870093901871995667,
            338064578630170283822037212269605496043572650\
                6343143038572*i + 33547126451073866168516\
                23808842170336793174147440300290760,
            1
        ],
        [
            305911177965252789833383994154493036996770951\
                0172993886444*i + 92565742612037761049150\
                4677137923505124932074435358013568,
            298080905845429690262707890078011680618909522\
                8301685298308*i + 34322372318479954444759\
                97009861549810901885286817560821549,
            995462354109678197105807353386077887495415751\
                245688653545*i + 192754851293894114748183\
                7522784332496852136496595867631368,
            1
        ],
        [
            218719111498740759726484473730219221250413037\
                4335649016672*i + 38471505181864174412729\
                40069955814813932425908079263797264,
            188149548951853215627001374103726727356196382\
                1880852849260*i + 85012682885534394906916\
                4281311480191554348089643106606951,
            208685434510057657110752947964880646853787126\
                6292849054213*i + 12020497251458084451890\
                30869106846459271749360294017240142,
            1
        ]
    ],
    [
        [
            146477681418900001186386142688798471674506225\
                2062842610827*i + 35247870364703646520982\
                48006865635670980904037259014207311,
            303107096550488297718593685076647376090184118\
                0865151334734*i + 13127078998957510162106\
                560481129931830554719058083966906,
            134058427310815187415911660344603688563489962\
                3831694175107*i + 14168154170108475631894\
                87197164230512137563777173228536319,
            1
        ],
        [
            117122087574744009448759036319539076483719329\
                6817013826379*i + 39454405214685111066247\
                7865457600116929217649321699053575,
            514148830212342255720693841275757820409275704\
                917360080435*i + 348554099774466102168288\
                3371128605669864415130509736034089,
            356862275263027428743359531475070698293929790\
                9171827009862*i + 16446149641929404181025\
                12034306237089908468705248427146529,
            1
        ],
        [
            119543528395684542902978787225554645401214762\
                9220540416182*i + 29909805247168385199087\
                41487732592391270140665192245111521,
            102651845758685297666442665487666256765979094\
                2087189146026*i + 28899457099247379430433\
                5318936269696248452518294207217636,
            201054976958576997259122402710304206150144501\
                3431697504787*i + 42552088611985916342332\
                55290260961792032991075601559394586,
            1
        ],
        [
            327600866484795251058117002246626731369086656\
                1590705122837*i + 39328190754315684625055\
                13238170139608307333301310535607064,
            356673059315439228808868109375830814345688763\
                2770392500629*i + 18828540995818559351941\
                31691041417935843653513207098391549,
            240542045254670176213964933831819112553058556\
                8565234938558*i + 14487580884365188221290\
                38573266690221326946786792342384705,
            1
        ]
    ],
    [
        [
            305911177965252789833383994154493036996770951\
                0172993886444*i + 92565742612037761049150\
                4677137923505124932074435358013568,
            298080905845429690262707890078011680618909522\
                8301685298308*i + 34322372318479954444759\
                97009861549810901885286817560821549,
            995462354109678197105807353386077887495415751\
                245688653545*i + 192754851293894114748183\
                7522784332496852136496595867631368,
            1
        ],
        [
            322976807103332165146148181113248196660064932\
                5392559718602*i + 32128018232436803164244\
                60650119693628695876941007409008035,
            129506428974323366430162360304434146718789665\
                1732402736102*i + 93769068037729199848950\
                1472940615413224522039795090905695,
            228497068144326080559926525740325440139944192\
                2369921834090*i + 23981631512212292618880\
                92150759916000269326600272170787172,
            1
        ],
        [
            305911177965252789833383994154493036996770951\
                0172993886444*i + 92565742612037761049150\
                4677137923505124932074435358013568,
            298080905845429690262707890078011680618909522\
                8301685298308*i + 34322372318479954444759\
                97009861549810901885286817560821549,
            995462354109678197105807353386077887495415751\
                245688653545*i + 192754851293894114748183\
                7522784332496852136496595867631368,
            1
        ],
        [
            981794237711585705231193519319833828701640042\
                099256443710*i + 306282757545837273046796\
                7476809436515599659793145226019091,
            172792165701223891678072905354130085420771806\
                1670845508880*i + 17614844124879692803094\
                68367108718023235810745247500720365,
            190934883927252048530893619269279649807041839\
                5431345345895*i + 27950487162292482694080\
                39901498031266972853105705182345932,
            1
        ]
    ],
    [
        [
            318061222727715218215094343579627702176550032\
                89375539139*i + 1020719658316342498347286\
                503110742505967417215854384643316,
            313011583926264772514846877546187774690030488\
                3925282354478*i + 50511150760349870752306\
                3753765033348657300943396885453085,
            422057279207623869511260778469104421613550943\
                2729031297689*i + 21052163327757514892731\
                5881629447112886313269875222127189,
            1
        ],
        [
            254189867054532339101846669674285727138806408\
                6970879515442*i + 10583126258379777096348\
                38840884875080817693628368747427569,
            218759986458508666349438519622066986324977456\
                810349181996*i + 333545735753950379126666\
                2593941254062416870093901871995667,
            338064578630170283822037212269605496043572650\
                6343143038572*i + 33547126451073866168516\
                23808842170336793174147440300290760,
            1
        ],
        [
            305911177965252789833383994154493036996770951\
                0172993886444*i + 92565742612037761049150\
                4677137923505124932074435358013568,
            298080905845429690262707890078011680618909522\
                8301685298308*i + 34322372318479954444759\
                97009861549810901885286817560821549,
            995462354109678197105807353386077887495415751\
                245688653545*i + 192754851293894114748183\
                7522784332496852136496595867631368,
            1
        ],
        [
            218719111498740759726484473730219221250413037\
                4335649016672*i + 38471505181864174412729\
                40069955814813932425908079263797264,
            188149548951853215627001374103726727356196382\
                1880852849260*i + 85012682885534394906916\
                4281311480191554348089643106606951,
            208685434510057657110752947964880646853787126\
                6292849054213*i + 12020497251458084451890\
                30869106846459271749360294017240142,
            1
        ]
    ],
    [
        [
            184948743102040794029302914922554579178770815\
                2625022733953*i + 29858993325538022391130\
                27539872569697527974789100881206249,
            285865391392243663258409488846957215066374028\
                7136499090572*i + 37746203869708977707505\
                09901981425855794861834486163917180,
            206006293821273059842121329750801927082928763\
                4936681633561*i + 26296489466195147939758\
                95848475500107151871259450546248699,
            1
        ],
        [
            286056598381284357757665429458887194941610197\
                1161222110243*i + 40137293811872286425602\
                34039716628666469097630235500710922,
            280308989661314215851641807443214002733817717\
                2044013844258*i + 49497548557853645302710\
                9687989363439776767819335928696553,
            113949281461202904760802438435732988004503512\
                4028650189315*i + 27879442209869141900589\
                89143864295372213551193413001848030,
            1
        ],
        [
            307733445433194636320426173972752251265791405\
                6747519142446*i + 34821075829171034572685\
                22599717880952429791658451506699094,
            152087222098031690089873254097552908282844690\
                9311176322364*i + 32851895318940802638455\
                71007101660966802759862687750849666,
            301129529489906756387432199900257297165639550\
                752782883276*i + 243826125444813489694850\
                3347414530453688240456450512850521,
            1
        ],
        [
            309196022323656640639748037559991660592319953\
                5310583580183*i + 19189936251645161680759\
                16935048687176202942956543232017079,
            335717806174151690331052604533285563739156070\
                956684147336*i + 201427530431048556761037\
                7847715320183759000952840304520672,
            115517796480228325758486343815842971011373404\
                0523304017775*i + 13915436305775393943878\
                49381982051751068370153675360350493,
            1
        ]
    ]
];
