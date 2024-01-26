//clear;

//////////////////////////////////
// Parameters for secp = 128
//////////////////////////////////

e:=75;
p:=5*2^4*3^e-1;
Fp:=GF(p);
Fp2<i>:=ExtensionField<Fp,x|x^2+1>;
_<x>:=PolynomialRing(Fp2);

//Top Strategy: [ 42, 33 ]
//217 TripleKummer and 255 Isogeny33Eval ==  16406 total units
strategy:=[ 0, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 11, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 18, 19, 20, 21, 21, 21, 21, 21, 22, 23, 24, 25, 26, 27, 27, 28, 29, 30, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33 ];


//load "setup.m";
//K,consts:=setup(f,4,e);

K:=[
    [
        31686932248522568585099010551058246465*i +
            35520209765133275254634260044078574569,
        20331847229163256069377052247152064943*i +
            3306988340664736660887945659699987863,
        22586265731572579411462439905328603410*i +
            41030671852390707390145502087579990483,
        22236804450089770444608186179282207377*i +
            26194855615263923747845216317273305162
    ],
    [
        3098710989611319931233055364401186926*i +
            8444121715157197025460480040936705785,
        29999998035773188660322103762422128620*i +
            32384036501083423447454466561420274176,
        26360190495763027781951331648627791275*i +
            38681085635918623588927041545551618844,
        1
    ],
    [
        3125994806732072787997797909573767602*i +
            3727608937064124124347276817457540905,
        27123413076082936159906870415209176647*i +
            3822346034488668488265751226425806907,
        11559178644001289007523643502461966777*i +
            25887153711927930733735842280933304298,
        9573805690821319016979766594414593065*i +
            36677350664177300915233678144733447326
    ],
    [
        20355445429154945449818840843682710260*i +
            45096551795483903327109610645241028403,
        34790360785908048362336283016694159572*i +
            40625599886475017540943283583583676988,
        20418640220495032486644558574555145811*i +
            3968215743907707110495864410374407492,
        6216465592971928683452921609852818200*i +
            13262395051421026286921225941348402289
    ]
];

gens:=[
    [
        [
            43747157381526401710080290936648948612*i +
                26241717368541231167805856286859667276,
            6561135426902168995707978978450755745*i +
                38174863830396770118023824115912973193,
            13083296322870436186946623529863302654*i +
                36420878148855316611746219478328071779,
            1
        ],
        [
            26017848369273854891534564300429894827*i +
                36676818428536422344233566583202052022,
            28993018911095840603081894886567061672*i +
                17673942823181231305690266927363793938,
            8085086699375893122321945787696078887*i +
                42168835235905053170324077551001091713,
            1
        ],
        [
            43747157381526401710080290936648948612*i +
                26241717368541231167805856286859667276,
            6561135426902168995707978978450755745*i +
                38174863830396770118023824115912973193,
            13083296322870436186946623529863302654*i +
                36420878148855316611746219478328071779,
            1
        ],
        [
            3842532230009761877442210281735915931*i +
                47249715520067245778288080852723646794,
            27111426523207847252382198305655395120*i +
                32887728851243603485101244367333217676,
            30147596569302746221026570453015264389*i +
                42400348373172067630466049051193399804,
            1
        ]
    ],
    [
        [
            22362523227295494999487613254003371191*i +
                40324036105970123323865554966269640783,
            19097807297855281054843266975939185485*i +
                36920903781299606407277564548234895954,
            29113494440851253582105957068553116944*i +
                20474262318025172976017838934192316634,
            1
        ],
        [
            26870744238559105342132957197406495709*i +
                19998193979596866026292131503514428424,
            30141947869442169880436854389572896440*i +
                40154269880464781424633114172934444344,
            41823262077358134478750394186982133733*i +
                41895260031237614229006255260418284212,
            1
        ],
        [
            43747157381526401710080290936648948612*i +
                26241717368541231167805856286859667276,
            6561135426902168995707978978450755745*i +
                38174863830396770118023824115912973193,
            13083296322870436186946623529863302654*i +
                36420878148855316611746219478328071779,
            1
        ],
        [
            19862228011425530645938514302151448461*i +
                3638489084140708302869569620581329901,
            26026975246392674674673048288965569571*i +
                27741450623886400516348391167634783481,
            44573389521525192055180000909231628416*i +
                41299206594704359234786361275550404530,
            1
        ]
    ],
    [
        [
            9299738932542353338788903268498215441*i +
                27956851874165671667079999146930684914,
            16830568136788711975715135717388963307*i +
                27364667912166024359752791157929144521,
            38563834563589741937609295599229369650*i +
                415114702101804947016660138194379464,
            1
        ],
        [
            43075515735002863895750856201320654402*i +
                44842973665826720204050911871946037872,
            48127896063334150900720020402989950222*i +
                22980761225138136490472104424749135071,
            18750804503107512071116912653261147139*i +
                48502986806082208230227580378621634646,
            1
        ],
        [
            9113146088277712276646924802052797323*i +
                33681872308926658638483635106779220144,
            499400009144705979755811831967409691*i +
                23983709393797865637103575835004260990,
            27020190739278148048324249485453325383*i +
                47496133113531672895729379330584820099,
            1
        ],
        [
            38483494412534343298348956490073431959*i +
                33555590153458334073131057773926404705,
            15481202737615192664609837620256821618*i +
                25343910588037069769551763119916790045,
            81101214784840708460655125284269818*i +
                22216438390931196498949677552471430486,
            1
        ]
    ],
    [
        [
            43747157381526401710080290936648948612*i +
                26241717368541231167805856286859667276,
            6561135426902168995707978978450755745*i +
                38174863830396770118023824115912973193,
            13083296322870436186946623529863302654*i +
                36420878148855316611746219478328071779,
            1
        ],
        [
            26017848369273854891534564300429894827*i +
                36676818428536422344233566583202052022,
            28993018911095840603081894886567061672*i +
                17673942823181231305690266927363793938,
            8085086699375893122321945787696078887*i +
                42168835235905053170324077551001091713,
            1
        ],
        [
            43747157381526401710080290936648948612*i +
                26241717368541231167805856286859667276,
            6561135426902168995707978978450755745*i +
                38174863830396770118023824115912973193,
            13083296322870436186946623529863302654*i +
                36420878148855316611746219478328071779,
            1
        ],
        [
            21779184616727075324068541941612981108*i +
                31551044688319565221916094077536904921,
            46941988125052022297776532058123313433*i +
                15418922421218175479701741144137174261,
            26471465297453118459706227409930981708*i +
                47282439931860442100294431739326232259,
            1
        ]
    ],
    [
        [
            22362523227295494999487613254003371191*i +
                40324036105970123323865554966269640783,
            19097807297855281054843266975939185485*i +
                36920903781299606407277564548234895954,
            29113494440851253582105957068553116944*i +
                20474262318025172976017838934192316634,
            1
        ],
        [
            26870744238559105342132957197406495709*i +
                19998193979596866026292131503514428424,
            30141947869442169880436854389572896440*i +
                40154269880464781424633114172934444344,
            41823262077358134478750394186982133733*i +
                41895260031237614229006255260418284212,
            1
        ],
        [
            43747157381526401710080290936648948612*i +
                26241717368541231167805856286859667276,
            6561135426902168995707978978450755745*i +
                38174863830396770118023824115912973193,
            13083296322870436186946623529863302654*i +
                36420878148855316611746219478328071779,
            1
        ],
        [
            19862228011425530645938514302151448461*i +
                3638489084140708302869569620581329901,
            26026975246392674674673048288965569571*i +
                27741450623886400516348391167634783481,
            44573389521525192055180000909231628416*i +
                41299206594704359234786361275550404530,
            1
        ]
    ],
    [
        [
            30356713000005612104935196081413378723*i +
                17542472505181092868334850177502789757,
            48653043090066322640143666732190430874*i +
                38378486272982754020204225830755803848,
            3703745984028737183614272301420426595*i +
                2658051042618291599552900801448028488,
            1
        ],
        [
            32582447945321890609541351757967323304*i +
                10861373468359740100474335408141327740,
            36699879744921546007246536007159593829*i +
                48289509292958082487695371603241493752,
            36268914632699322986150320586512316416*i +
                32742407928187178786142235474921515350,
            1
        ],
        [
            9551389466914005203282214691842350763*i +
                42778156791834635626221097329508057843,
            17740819460993838887635184077163957791*i +
                37402607879379951006698065677671139401,
            33185057177805817166244480758217617*i +
                34639045187697442882229511719566726281,
            1
        ],
        [
            2666156295530671340435701844507070772*i +
                26479509992039125835011500002766314023,
            32398811296039949696975440197320164054*i +
                11037391014527724186778187694631384751,
            10013302137492149842497990041835234651*i +
                4428308325185875008447085523942770545,
            1
        ]
    ]
];
