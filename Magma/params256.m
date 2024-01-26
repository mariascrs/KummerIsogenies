//clear;

//////////////////////////////////
// Parameters for secp = 256
//////////////////////////////////

e:=154;
p:=11*2^4*3^e-1;
Fp:=GF(p);
Fp2<i>:=ExtensionField<Fp,x|x^2+1>;
_<x>:=PolynomialRing(Fp2);

//Top Strategy for A: [ 86, 68 ]
//517 TripleKummer and 613 Isogeny33Eval ==  39262 total units
strategy:=[ 0, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 11, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 18, 19, 20, 21, 21, 21, 21, 21, 22, 23, 24, 25, 26, 27, 27, 28, 29, 30,
31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 35, 36, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 48, 48, 48, 48, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 58, 59, 60, 61, 62, 63, 63, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, 67, 68 ];

//load "setup.m";
//K,consts:=setup(f,4,e);

K:=[
    [
        4678573165223611166191414910975215918572812371594\
            782534368734734098854043250*i +
            355986250236857836555252230457179860659717713\
            8986167221301954321745296808594,
        1168891729736991904617964593612110349543677529368\
            332469956309163159443165679*i +
            349603655154783071947300320950988884800554123\
            6148240943572772597616458020749,
        3513150309085340847966775316241562368993723653012\
            959871656377420155937837368*i +
            142897537876880241435918339082414891231100180\
            3789399320273659753500928827733,
        7655498206593612554094565346136878213774476438320\
            24375532439103727834159703*i +
            302663991579393343299400685063454203756403489\
            462106677032766179672939429920
    ],
    [
        4728440471006727476411980691449015196753945522963\
            513540097170067028163424282*i +
            297568284756131312275382349712556497683141461\
            2017703525326344552131550004465,
        1912275635628200946118636771470748023768337058487\
            598556154128441324436885624*i +
            409943634060718811750565664263694979615588173\
            441232627405268296247484493632,
        2939637034561908365092371380929963557449521471002\
            961848496246197232733133555*i +
            258564470390762888992069129640057932969558678\
            1197598242856848560500823021434,
        1
    ],
    [
        2989240182966962120793451580150395602062346894173\
            993381675540181681029450943*i +
            821123245637948400503817026086539554489731710\
            251796305306449503789678805534,
        5246563335255365739716205566875689057500757202373\
            31119355505011011876013372*i +
            147297190414073751399666984689345584107862576\
            3016321106759218440474909316628,
        9347573150651304729324203469859806621040628712785\
            97709061877624325296611861*i +
            285126390307695274416813519046517816016418408\
            6074907860054083788540451795458,
        2193998835650933377115791302795369073395990228674\
            377053191523968912644008378*i +
            120206249290688819065541161336553566283936288\
            4664345696942601335746648638777
    ],
    [
        2271198461216529913612535750462965070479073760585\
            9433632638949875846134909*i +
            474213088840638336577733649534737867372034710\
            1930370416738393922825891016053,
        1671014632486752538846854387418584090400878897296\
            125417343362485246338232749*i +
            262277545319529938800171074979944455530016303\
            6763271797764414028308186382564,
        7258673201530959973273509470130499510789166166946\
            67937386086874634125372095*i +
            614302281346236418203994586945310610853737414\
            599070233157243342482936503055,
        4102466614870542113706193000953856113731337093221\
            477922613881963879882499832*i +
            139077965556171333766773539051942360111005654\
            1828448634899978397417656858240
    ]
];

gens:=[
    [
        [
            137616500475656160374077082987018956809149343\
                1017855879991415081237147058829*i +
                33834032349533466972486122183771725667545\
                72898737116734748175530796716929195,
            525119705218659805031616495819801017999080734\
                0215611953107899381033546554493*i +
                42710709820204893285158568421638462958025\
                30994948485691300078642078762706881,
            420513488713173864344615105321563047858609387\
                4592757369223417592173774611338*i +
                41939885523199301600085345317072562526782\
                4980734585190272686533030262175528,
            1
        ],
        [
            579542042162802426171269337201261848198430284\
                505870725312883170834233769078*i +
                32841998746581188794642069577173485240329\
                81005043278841254843911624791009951,
            126846261887949014635387469884953403086879159\
                8684800480090618041561112746205*i +
                34759771930096054098375991842058438815865\
                2467259259274325419871051338791193,
            173464330882782993215771909723623646647856650\
                7465493239720542272551214058298*i +
                45439049321741350233457449274050835364175\
                58223542090707528435738097589561127,
            1
        ],
        [
            137616500475656160374077082987018956809149343\
                1017855879991415081237147058829*i +
                33834032349533466972486122183771725667545\
                72898737116734748175530796716929195,
            525119705218659805031616495819801017999080734\
                0215611953107899381033546554493*i +
                42710709820204893285158568421638462958025\
                30994948485691300078642078762706881,
            420513488713173864344615105321563047858609387\
                4592757369223417592173774611338*i +
                41939885523199301600085345317072562526782\
                4980734585190272686533030262175528,
            1
        ],
        [
            333758503155534589088283575578825046205330852\
                0185763832720645260063463321490*i +
                15906869566635889708834398925864007004882\
                84468840514339519425784963007511735,
            273541575712586703601717909310105117405958994\
                8421079092308272756924428551932*i +
                47917904395693981532226197316743029757222\
                05181758689769433696373387734833668,
            102597683196469815708692367909988743226411044\
                0297211457004193469021795663639*i +
                77526486871544609446086159918010496486301\
                1775139658514671052071321997753472,
            1
        ]
    ],
    [
        [
            240269971580752781787620761559317400157865356\
                5623526101413771605423147458104*i +
                16028431277329029968402678198720196778586\
                07953874330820208017550442083336186,
            460131421545352350585144862264752584722791422\
                9723196503353108651300525741351*i +
                39560601992073345196559881865144434833756\
                60071806289996665401210877462956563,
            555481850569623666570836305577974716567787949\
                111205935056178210785825659402*i +
                32817944575265054606001867579497469509532\
                86356236360992580999443489652601000,
            1
        ],
        [
            105912155810551560006762391756269537101576745\
                1509100731750243190740436586427*i +
                21826145009647313037172807403164284191259\
                66574296894268487544228736313871123,
            240909203370213136537290019505930781527871225\
                2711444197689600554384687435245*i +
                21593543804422344184768109530613888583283\
                02521295385984471529209326945839887,
            443395193159796856978324698140866644528709195\
                0279615200708771823852134073350*i +
                45819374217630102863779617225909438520490\
                23291609391870600155898276688436986,
            1
        ],
        [
            137616500475656160374077082987018956809149343\
                1017855879991415081237147058829*i +
                33834032349533466972486122183771725667545\
                72898737116734748175530796716929195,
            525119705218659805031616495819801017999080734\
                0215611953107899381033546554493*i +
                42710709820204893285158568421638462958025\
                30994948485691300078642078762706881,
            420513488713173864344615105321563047858609387\
                4592757369223417592173774611338*i +
                41939885523199301600085345317072562526782\
                4980734585190272686533030262175528,
            1
        ],
        [
            189019504931938669070474907706577215811047362\
                8633495112842003232998554544484*i +
                19402555260888319348471877480364445195216\
                51626661847144538641517165378127009,
            501823114080860090686128237894879082606936802\
                0236578373110242386447549163111*i +
                52175994740258343155287417431154778058461\
                70575096075382400659615380389506421,
            482880093991762254685305430485857741131400263\
                842069026364636666153524561441*i +
                43964815951740189505902364164915554811295\
                90749506565718137222745384656355412,
            1
        ]
    ],
    [
        [
            281981603489865091550918561613785117614741967\
                1826355395650250739575849196657*i +
                28609739638085146763568869832062607669510\
                04068529889260606030645674519286995,
            464326603390220084953781251837448092247438809\
                1842289586199526712445755109713*i +
                85763908900607251060094939865960265903243\
                2153552617246651154979745302684975,
            481832788103125215631138874445169440321086745\
                5168505165565638189324669502818*i +
                12984422504986231812996740343344319305419\
                42291598565696658733325423110100764,
            1
        ],
        [
            339898402586771283955431191507492992840402383\
                4496605795001738389513015728588*i +
                20834105446731375292064434506574863979626\
                24375820018687519442985530568437981,
            433005113009433846689277945607380387757901756\
                239692361105500804012349775648*i +
                99403939298988241928206568312529513283948\
                8870566584261301581965927164939833,
            405543134801973686709588274016103024222400686\
                1768520103326339821071442766238*i +
                31137172593702220939741744228299327741996\
                0182785320306536029993696734601507,
            1
        ],
        [
            211594254779076511306521425393201237649218009\
                2928019572271680600340396795179*i +
                11190221530138884927921382184967359541742\
                01421587948038840506597127496767836,
            505434423112799658524915892006053201029106876\
                9294788372743279677362999124312*i +
                14456480573955021963795636172194020993081\
                64987691734451755542621933928420258,
            401883661784699731836805311562558807086942059\
                2769063788009194806178159953524*i +
                34856434032712210562939354168533963484528\
                4047135877886828048955027578933034,
            1
        ],
        [
            409288924591831161972607316600856054336095025\
                97038279653683965625827359266*i +
                21045279501825270823428689951756096483942\
                34823644155791614641317025609275521,
            218611599867909298348967009455239655331405417\
                9588945149336245984577246120072*i +
                38013525975138517760141132880113539563201\
                69013427895964965707418990507769584,
            110798369352711866768053327743959550129755648\
                0226272512060134885496858265959*i +
                18306842680247535631951852490875138981712\
                40419829238072713051288408606459310,
            1
        ]
    ],
    [
        [
            137616500475656160374077082987018956809149343\
                1017855879991415081237147058829*i +
                33834032349533466972486122183771725667545\
                72898737116734748175530796716929195,
            525119705218659805031616495819801017999080734\
                0215611953107899381033546554493*i +
                42710709820204893285158568421638462958025\
                30994948485691300078642078762706881,
            420513488713173864344615105321563047858609387\
                4592757369223417592173774611338*i +
                41939885523199301600085345317072562526782\
                4980734585190272686533030262175528,
            1
        ],
        [
            579542042162802426171269337201261848198430284\
                505870725312883170834233769078*i +
                32841998746581188794642069577173485240329\
                81005043278841254843911624791009951,
            126846261887949014635387469884953403086879159\
                8684800480090618041561112746205*i +
                34759771930096054098375991842058438815865\
                2467259259274325419871051338791193,
            173464330882782993215771909723623646647856650\
                7465493239720542272551214058298*i +
                45439049321741350233457449274050835364175\
                58223542090707528435738097589561127,
            1
        ],
        [
            137616500475656160374077082987018956809149343\
                1017855879991415081237147058829*i +
                33834032349533466972486122183771725667545\
                72898737116734748175530796716929195,
            525119705218659805031616495819801017999080734\
                0215611953107899381033546554493*i +
                42710709820204893285158568421638462958025\
                30994948485691300078642078762706881,
            420513488713173864344615105321563047858609387\
                4592757369223417592173774611338*i +
                41939885523199301600085345317072562526782\
                4980734585190272686533030262175528,
            1
        ],
        [
            467148894121116532588867684687350335396478902\
                6861377532337993564975260800509*i +
                19952543825065630647875625086329103440133\
                44251717825770894843800195199675419,
            238032624308363821887935641166869508279106923\
                0583555719054110709077598687804*i +
                23188111026261353858966616075780633602271\
                33086353714145073661437948897725334,
            523035147785328939339898708249961991322805567\
                2118597351381575498273901360502*i +
                10650634736894115418781422968699228882932\
                69669150003549371988551140084938645,
            1
        ]
    ],
    [
        [
            240269971580752781787620761559317400157865356\
                5623526101413771605423147458104*i +
                16028431277329029968402678198720196778586\
                07953874330820208017550442083336186,
            460131421545352350585144862264752584722791422\
                9723196503353108651300525741351*i +
                39560601992073345196559881865144434833756\
                60071806289996665401210877462956563,
            555481850569623666570836305577974716567787949\
                111205935056178210785825659402*i +
                32817944575265054606001867579497469509532\
                86356236360992580999443489652601000,
            1
        ],
        [
            105912155810551560006762391756269537101576745\
                1509100731750243190740436586427*i +
                21826145009647313037172807403164284191259\
                66574296894268487544228736313871123,
            240909203370213136537290019505930781527871225\
                2711444197689600554384687435245*i +
                21593543804422344184768109530613888583283\
                02521295385984471529209326945839887,
            443395193159796856978324698140866644528709195\
                0279615200708771823852134073350*i +
                45819374217630102863779617225909438520490\
                23291609391870600155898276688436986,
            1
        ],
        [
            137616500475656160374077082987018956809149343\
                1017855879991415081237147058829*i +
                33834032349533466972486122183771725667545\
                72898737116734748175530796716929195,
            525119705218659805031616495819801017999080734\
                0215611953107899381033546554493*i +
                42710709820204893285158568421638462958025\
                30994948485691300078642078762706881,
            420513488713173864344615105321563047858609387\
                4592757369223417592173774611338*i +
                41939885523199301600085345317072562526782\
                4980734585190272686533030262175528,
            1
        ],
        [
            189019504931938669070474907706577215811047362\
                8633495112842003232998554544484*i +
                19402555260888319348471877480364445195216\
                51626661847144538641517165378127009,
            501823114080860090686128237894879082606936802\
                0236578373110242386447549163111*i +
                52175994740258343155287417431154778058461\
                70575096075382400659615380389506421,
            482880093991762254685305430485857741131400263\
                842069026364636666153524561441*i +
                43964815951740189505902364164915554811295\
                90749506565718137222745384656355412,
            1
        ]
    ],
    [
        [
            230944989805093279076426467728719051591235378\
                6853248695671772220723123402094*i +
                43223104020935146497878626303126091514920\
                64952335134707916247956836815215063,
            377308630446792973160731130713020110909969267\
                0398324283887554291582173369074*i +
                51407086518467014009617567493045659681462\
                95226194705259632722048900142188697,
            390521592235882817704835037694522447599557097\
                8785592996281267575104529084348*i +
                22417599717563479293491441728850366492224\
                22831716507912959122997209579060142,
            1
        ],
        [
            334512118220745196733139263542884544434706758\
                6233583785786928897336695943484*i +
                63941638188556164319427584759244756890447\
                0606945588159729825465433530645133,
            326554418352137482065390308843393837887413036\
                4548929748967048903196423698972*i +
                61082364613604534151442019979429546212802\
                3990735767428600954225632017189627,
            504736504941145334511472861362343557261491967\
                6737666950839763735281119599202*i +
                17500762590660105278540056266661836271035\
                76274269241800985300235686618245077,
            1
        ],
        [
            371262110831902831827311006018545519668443551\
                4886380303927987859878667469200*i +
                49813623676441754550831329606626484105376\
                84627576946941352672777094249425464,
            399196136598168938231484142052712194345746898\
                3786578209334865047303366409884*i +
                28858945160541648606148672770340803566729\
                77998096241430694503132557849199540,
            368681719200054431852317033197244551255610250\
                515080721461787720400455031871*i +
                17492037278257713698481752099389400726847\
                78238794837950740513883183436798556,
            1
        ],
        [
            226241293011868742547653133697286229915120494\
                3468042982844910671800330448046*i +
                48837133418599447394535213638564405280916\
                84453420134674127664991799296069901,
            273024302071214042704998119834361675636872457\
                0624712653971619100624401293294*i +
                42404774241216380696057612081293369660254\
                80392803335717671528719795463162521,
            265732380548433307988157147919509302478914821\
                0843568011448273757499564895469*i +
                32907889833194674777085675008712220710888\
                51041988124210312304201079035206763,
            1
        ]
    ]
];
