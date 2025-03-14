$offlisting
$offdigit

EQUATIONS
	c1
	c2
	c3
	c4
	c5
	c6
	c7
	c8
	c9
	c10
	c11
	c12
	c13
	c14
	c15
	c16_hi
	c17_lo
	c18_hi
	c19_lo
	c20_hi
	c21_hi
	c22_hi
	c23_hi
	c24_hi
	c25_hi
	c26_hi
	c27_hi
	c28
	c29_hi
	c30
	c31_hi
	c32_hi
	c33_hi
	c34_hi
	c35_hi
	c36_hi
	c37_hi
	c38_hi
	c39_hi
	c40_hi
	c41_hi
	c42_hi
	c43_hi
	c44_hi
	c45_hi
	c46_hi
	c47_hi
	c48_hi
	c49_hi
	c50_hi
	c51_hi
	c52_lo
	c53_lo
	c54_lo
	c55_lo
	c56_lo
	c57_lo
	c58_lo
	c59_lo
	c60_lo
	c61_lo
	c62_lo
	c63_lo
	c64_lo
	c65_lo
	c66_lo
	c67_lo
	c68_lo
	c69_lo
	c70_lo
	c71_lo
	c72_lo
	c73_lo
	c74_lo
	c75_lo
	c76_lo
	c77_lo
	c78_lo
	c79_lo
	c80_lo
	c81_lo
	c82_lo
	c83_lo
	c84_hi
	c85_hi
	c86_hi
	c87_hi
	c88_hi
	c89_hi
	c90_hi
	c91_hi
	c92_hi
	c93_hi
	c94_hi
	c95_hi
	c96_hi
	c97_hi
	c98_hi
	c99_hi
	c100_hi
	c101_hi
	c102_hi
	c103_hi
	c104_lo
	c105_lo
	c106_lo
	c107_lo
	c108_lo
	c109_lo
	c110_lo
	c111_lo
	c112_lo
	c113_lo
	c114_lo
	c115_lo
	c116_lo
	c117_lo
	c118_lo
	c119_lo
	c120_lo
	c121_lo
	c122_lo
	c123_lo
	c124_lo
	c125_lo
	c126_lo
	c127_lo
	c128_lo
	c129_lo
	c130_lo
	c131_lo
	c132_lo
	c133_lo
	c134_lo
	c135_lo
	c136;

BINARY VARIABLES
	x25
	x26
	x27
	x28
	x29
	x30
	x31
	x32;

POSITIVE VARIABLES
	x1
	x2
	x3
	x4
	x5
	x6
	x7
	x8
	x9
	x10
	x11
	x12
	x13
	x14
	x15
	x16
	x17
	x18
	x19
	x20
	x21
	x22
	x23
	x24;

VARIABLES
	GAMS_OBJECTIVE
	;


c1.. -x1 + exp(x2) =e= 1 ;
c2.. -x3 + exp(0.833333333333333*x4) =e= 1 ;
c3.. -x5 + 1.5*x6 + x7 =e= 0 ;
c4.. 1.25*x8 - x9 + 1.25*x10 =e= 0 ;
c5.. x11 + (-2)*x12 =e= 0 ;
c6.. -x13 + exp(0.66666666666666696*x14) =e= 1 ;
c7.. -x15 + exp(x16) =e= 1 ;
c8.. -x7 - x17 + exp(x18) =e= 1 ;
c9.. x9 - x13 - x15 =e= 0 ;
c10.. -x6 - x12 + x17 - x19 =e= 0 ;
c11.. x20 - x8 - x11 =e= 0 ;
c12.. x2 + x4 - x21 - x20 =e= 0 ;
c13.. x21 - x22 - x5 =e= 0 ;
c14.. -x14 - x16 + x23 =e= 0 ;
c15.. -x10 + x23 - x24 =e= 0 ;
c16_hi.. x7 + (-0.8)*x17 =l= 0 ;
c17_lo.. 0 =l= x7 + (-0.4)*x17 ;
c18_hi.. x8 + (-5)*x10 =l= 0 ;
c19_lo.. 0 =l= x8 + (-2)*x10 ;
c20_hi.. (-10)*x25 + x1 =l= 0 ;
c21_hi.. (-10)*x26 + x3 =l= 0 ;
c22_hi.. (-10)*x27 + x6 =l= 0 ;
c23_hi.. (-10)*x28 + x8 + x10 =l= 0 ;
c24_hi.. (-10)*x29 + x11 =l= 0 ;
c25_hi.. (-10)*x30 + x13 =l= 0 ;
c26_hi.. (-10)*x31 + x15 =l= 0 ;
c27_hi.. (-10)*x32 + x7 + x17 =l= 0 ;
c28.. x26 + x25 =e= 1 ;
c29_hi.. x28 + x29 =l= 1 ;
c30.. -x28 + x30 + x31 =e= 0 ;
c31_hi.. x27 - x32 =l= 0 ;
c32_hi.. 1.07217*x2 - x1 + 0.32149*x2*x2 =l= 0.0062857170715240329 ;
c33_hi.. 1.3354699999999999*x2 - x1 + 0.49213*x2*x2 =l= 0.16227611100764305 ;
c34_hi.. 16.39458*x2 - x1 + (-2.3521399999999999)*x2*x2 =l= 16.998306821459149 ;
c35_hi.. 0.92996*x2 - x1 + 0.25018*x2*x2 =l= 0.0054753055833179909 ;
c36_hi.. 1.6226499999999999*x2 - x1 + 0.78822*x2*x2 =l= 0.77224574866568019 ;
c37_hi.. 1.7855399999999999*x2 - x1 + 1.1486099999999999*x2*x2 =l= 1.8880853788009102 ;
c38_hi.. 1.21192*x2 - x1 + 0.40455*x2*x2 =l= 0.05902572054163768 ;
c39_hi.. 1.74783*x2 - x1 + 1.02457*x2*x2 =l= 1.4673109365414025 ;
c40_hi.. 1.8069299999999999*x2 - x1 + 1.2734399999999999*x2*x2 =l= 2.3445850522111122 ;
c41_hi.. 0.94223999999999997*x2 - x1 + 0.25589*x2*x2 =l= 0.003746519007379656 ;
c42_hi.. 1.32915*x2 - x1 + 0.48726*x2*x2 =l= 0.15541354963579557 ;
c43_hi.. 1.54949*x2 - x1 + 0.69394999999999996*x2*x2 =l= 0.54151527996159676 ;
c44_hi.. 0.83681*x2 - x1 + 0.20935*x2*x2 =l= 0.028454819305588797 ;
c45_hi.. 0.18165999999999999*x2 - x1 + (-0.0030599999999999998)*x2*x2 =l= 0.49988143537200946 ;
c46_hi.. 1.4432*x2 - x1 + 0.58350999999999997*x2*x2 =l= 0.31235872913507889 ;
c47_hi.. 1.15855*x2 - x1 + 0.37108*x2*x2 =l= 0.03188572381839749 ;
c48_hi.. 0.47749999999999998*x2 - x1 + 0.08142*x2*x2 =l= 0.24453664148147936 ;
c49_hi.. 0.67013*x2 - x1 + 0.14605*x2*x2 =l= 0.10861676964144529 ;
c50_hi.. 1.6913899999999999*x2 - x1 + 0.90002*x2*x2 =l= 1.0819895566225832 ;
c51_hi.. 1.45363*x2 - x1 + 0.59333*x2*x2 =l= 0.33072808007793064 ;
c52_lo.. -0.9904973696030137 =l= 1.8118099999999999*x2 - x1 + 0.44374*x2*x2 ;
c53_lo.. -0.02506234199226598 =l= 0.73728*x2 - x1 + 1.22237*x2*x2 ;
c54_lo.. -0.7201633481854199 =l= (-0.82521)*x2 - x1 + 1.82983*x2*x2 ;
c55_lo.. -0.5357730551333306 =l= 1.7404299999999999*x2 - x1 + 0.59311*x2*x2 ;
c56_lo.. -0.81965717127234328 =l= 1.79794*x2 - x1 + 0.49337999999999999*x2*x2 ;
c57_lo.. -0.035857144971769728 =l= 1.2666999999999999*x2 - x1 + 0.95494999999999997*x2*x2 ;
c58_lo.. -0.15270269348652499 =l= 1.4932399999999999*x2 - x1 + 0.81247*x2*x2 ;
c59_lo.. -3.6055178888362125 =l= (-4.4785899999999996)*x2 - x1 + 2.93522*x2*x2 ;
c60_lo.. -0.42192210362016858 =l= 1.69659*x2 - x1 + 0.64349*x2*x2 ;
c61_lo.. -0.26822978488652116 =l= 1.60554*x2 - x1 + 0.72743999999999998*x2*x2 ;
c62_lo.. -0.5255100128153047 =l= 1.73708*x2 - x1 + 0.59735*x2*x2 ;
c63_lo.. -0.0008971398421537824 =l= 1.04518*x2 - x1 + 1.0744499999999999*x2*x2 ;
c64_lo.. -0.34000683628396167 =l= 1.65385*x2 - x1 + 0.6853399999999999*x2*x2 ;
c65_lo.. -0.9877356524845109 =l= 1.8117099999999999*x2 - x1 + 0.44447999999999999*x2*x2 ;
c66_lo.. -0.12598851513037945 =l= 1.4573*x2 - x1 + 0.83711999999999998*x2*x2 ;
c67_lo.. -0.088402210018736938 =l= 0.47697*x2 - x1 + 1.33668*x2*x2 ;
c68_lo.. -0.6381937896338756 =l= 1.76852*x2 - x1 + 0.5534599999999999*x2*x2 ;
c69_lo.. -0.41721716115337326 =l= 1.69442*x2 - x1 + 0.64575*x2*x2 ;
c70_lo.. -0.572037165389764 =l= 1.75143*x2 - x1 + 0.57854*x2*x2 ;
c71_lo.. -1.2563906103746745 =l= (-1.62568)*x2 - x1 + 2.09601*x2*x2 ;
c72_lo.. -1.29627273939581 =l= (-1.68114)*x2 - x1 + 2.11377*x2*x2 ;
c73_lo.. -0.9529371741829748 =l= 1.8101799999999999*x2 - x1 + 0.45394*x2*x2 ;
c74_lo.. -0.11449404815316866 =l= 1.4400299999999999*x2 - x1 + 0.84863*x2*x2 ;
c75_lo.. -0.137886983032669 =l= 1.47397*x2 - x1 + 0.82581*x2*x2 ;
c76_lo.. -0.0048935629851754037 =l= 0.88836*x2 - x1 + 1.1518699999999999*x2*x2 ;
c77_lo.. -0.25958721951594876 =l= 1.59885*x2 - x1 + 0.73294999999999999*x2*x2 ;
c78_lo.. -3.4242771018714926 =l= (-4.2767099999999996)*x2 - x1 + 2.87955*x2*x2 ;
c79_lo.. -1.0376086973113288 =l= (-1.9983599999999999)*x2 - x1 + 2.80798*x2*x2 ;
c80_lo.. -2.0713197054097789 =l= (-2.69211)*x2 - x1 + 2.4254899999999999*x2*x2 ;
c81_lo.. -2.1097540871859977 =l= (-2.7397)*x2 - x1 + 2.43968*x2*x2 ;
c82_lo.. -0.73732877714751477 =l= 1.78791*x2 - x1 + 0.51898*x2*x2 ;
c83_lo.. -0.25359709893064997 =l= 0.040189999999999997*x2 - x1 + 1.51377*x2*x2 ;
c84_hi.. 0.89347499999999969*x4 - x3 + 0.2232569444444443*x4*x4 =l= 0.0062857170715240329 ;
c85_hi.. 1.1128916666666662*x4 - x3 + 0.3417569444444442*x4*x4 =l= 0.16227611100764305 ;
c86_hi.. 13.662149999999997*x4 - x3 + (-1.6334305555555544)*x4*x4 =l= 16.998306821459149 ;
c87_hi.. 0.77496666666666636*x4 - x3 + 0.173736111111111*x4*x4 =l= 0.0054753055833179909 ;
c88_hi.. 1.3522083333333328*x4 - x3 + 0.5473749999999996*x4*x4 =l= 0.77224574866568019 ;
c89_hi.. 1.4879499999999994*x4 - x3 + 0.7976458333333327*x4*x4 =l= 1.8880853788009102 ;
c90_hi.. 1.0099333333333331*x4 - x3 + 0.28093749999999984*x4*x4 =l= 0.05902572054163768 ;
c91_hi.. 1.4565249999999994*x4 - x3 + 0.71150694444444396*x4*x4 =l= 1.4673109365414025 ;
c92_hi.. 1.5057749999999994*x4 - x3 + 0.8843333333333326*x4*x4 =l= 2.3445850522111122 ;
c93_hi.. 0.78519999999999968*x4 - x3 + 0.17770138888888878*x4*x4 =l= 0.003746519007379656 ;
c94_hi.. 1.1076249999999996*x4 - x3 + 0.3383749999999998*x4*x4 =l= 0.15541354963579557 ;
c95_hi.. 1.2912416666666662*x4 - x3 + 0.48190972222222189*x4*x4 =l= 0.54151527996159676 ;
c96_hi.. 0.6973416666666664*x4 - x3 + 0.14538194444444436*x4*x4 =l= 0.028454819305588797 ;
c97_hi.. 0.15138333333333326*x4 - x3 + (-0.0021249999999999984)*x4*x4 =l= 0.49988143537200946 ;
c98_hi.. 1.2026666666666663*x4 - x3 + 0.40521527777777749*x4*x4 =l= 0.31235872913507889 ;
c99_hi.. 0.96545833333333297*x4 - x3 + 0.25769444444444428*x4*x4 =l= 0.03188572381839749 ;
c100_hi.. 0.39791666666666653*x4 - x3 + 0.056541666666666636*x4*x4 =l= 0.24453664148147936 ;
c101_hi.. 0.5584416666666665*x4 - x3 + 0.10142361111111105*x4*x4 =l= 0.10861676964144529 ;
c102_hi.. 1.4094916666666661*x4 - x3 + 0.6250138888888885*x4*x4 =l= 1.0819895566225832 ;
c103_hi.. 1.2113583333333329*x4 - x3 + 0.41203472222222198*x4*x4 =l= 0.33072808007793064 ;
c104_lo.. -0.9904973696030137 =l= 1.509841666666666*x4 - x3 + 0.3081527777777776*x4*x4 ;
c105_lo.. -0.02506234199226598 =l= 0.6143999999999998*x4 - x3 + 0.848868055555555*x4*x4 ;
c106_lo.. -0.7201633481854199 =l= (-0.6876749999999997)*x4 - x3 + 1.270715277777777*x4*x4 ;
c107_lo.. -0.5357730551333306 =l= 1.4503583333333327*x4 - x3 + 0.4118819444444442*x4*x4 ;
c108_lo.. -0.81965717127234328 =l= 1.498283333333333*x4 - x3 + 0.34262499999999974*x4*x4 ;
c109_lo.. -0.035857144971769728 =l= 1.0555833333333329*x4 - x3 + 0.6631597222222217*x4*x4 ;
c110_lo.. -0.15270269348652499 =l= 1.2443666666666662*x4 - x3 + 0.56421527777777747*x4*x4 ;
c111_lo.. -3.6055178888362125 =l= (-3.7321583333333317)*x4 - x3 + 2.038347222222221*x4*x4 ;
c112_lo.. -0.42192210362016858 =l= 1.4138249999999994*x4 - x3 + 0.44686805555555525*x4*x4 ;
c113_lo.. -0.26822978488652116 =l= 1.3379499999999995*x4 - x3 + 0.5051666666666663*x4*x4 ;
c114_lo.. -0.5255100128153047 =l= 1.4475666666666662*x4 - x3 + 0.4148263888888886*x4*x4 ;
c115_lo.. -0.0008971398421537824 =l= 0.870983333333333*x4 - x3 + 0.7461458333333327*x4*x4 ;
c116_lo.. -0.34000683628396167 =l= 1.3782083333333328*x4 - x3 + 0.47593055555555519*x4*x4 ;
c117_lo.. -0.9877356524845109 =l= 1.5097583333333326*x4 - x3 + 0.3086666666666664*x4*x4 ;
c118_lo.. -0.12598851513037945 =l= 1.2144166666666663*x4 - x3 + 0.5813333333333329*x4*x4 ;
c119_lo.. -0.088402210018736938 =l= 0.39747499999999986*x4 - x3 + 0.92824999999999946*x4*x4 ;
c120_lo.. -0.6381937896338756 =l= 1.4737666666666662*x4 - x3 + 0.38434722222222195*x4*x4 ;
c121_lo.. -0.41721716115337326 =l= 1.4120166666666663*x4 - x3 + 0.4484374999999997*x4*x4 ;
c122_lo.. -0.572037165389764 =l= 1.4595249999999995*x4 - x3 + 0.40176388888888864*x4*x4 ;
c123_lo.. -1.2563906103746745 =l= (-1.3547333333333329)*x4 - x3 + 1.4555624999999992*x4*x4 ;
c124_lo.. -1.29627273939581 =l= (-1.4009499999999995)*x4 - x3 + 1.4678958333333325*x4*x4 ;
c125_lo.. -0.9529371741829748 =l= 1.5084833333333327*x4 - x3 + 0.31523611111111088*x4*x4 ;
c126_lo.. -0.11449404815316866 =l= 1.2000249999999995*x4 - x3 + 0.5893263888888884*x4*x4 ;
c127_lo.. -0.137886983032669 =l= 1.2283083333333329*x4 - x3 + 0.5734791666666663*x4*x4 ;
c128_lo.. -0.0048935629851754037 =l= 0.7402999999999997*x4 - x3 + 0.79990972222222168*x4*x4 ;
c129_lo.. -0.25958721951594876 =l= 1.3323749999999996*x4 - x3 + 0.50899305555555518*x4*x4 ;
c130_lo.. -3.4242771018714926 =l= (-3.5639249999999985)*x4 - x3 + 1.9996874999999987*x4*x4 ;
c131_lo.. -1.0376086973113288 =l= (-1.6652999999999993)*x4 - x3 + 1.9499861111111099*x4*x4 ;
c132_lo.. -2.0713197054097789 =l= (-2.2434249999999993)*x4 - x3 + 1.6843680555555545*x4*x4 ;
c133_lo.. -2.1097540871859977 =l= (-2.2830833333333325)*x4 - x3 + 1.6942222222222212*x4*x4 ;
c134_lo.. -0.73732877714751477 =l= 1.4899249999999995*x4 - x3 + 0.3604027777777775*x4*x4 ;
c135_lo.. -0.25359709893064997 =l= 0.033491666666666649*x4 - x3 + 1.051229166666666*x4*x4 ;
c136.. GAMS_OBJECTIVE =e= 8*x26 + 6*x27 + 10*x28 + 6*x29 + 7*x30 + 4*x31 + 5*x32 + 5*x25 + x1 + (-10)*x2 + x3 + (-15)*x4 + (-40)*x6 + 15*x7 + 15*x10 + 80*x17 + (-65)*x18 + 25*x13 + (-60)*x14 + 35*x15 + (-80)*x16 + (-35)*x19 + 122 ;

x2.up = 2;
x4.up = 2;
x6.up = 2;
x7.up = 1;
x10.up = 1;
x13.up = 2;
x14.up = 100;
x15.up = 2;
x16.up = 100;
x17.up = 2;
x18.up = 100;
x19.up = 3;

MODEL GAMS_MODEL /all/ ;
option solprint=off;
option limrow=0;
option limcol=0;
option solvelink=5;
SOLVE GAMS_MODEL USING minlp minimizing GAMS_OBJECTIVE;

Scalars MODELSTAT 'model status', SOLVESTAT 'solve status';
MODELSTAT = GAMS_MODEL.modelstat;
SOLVESTAT = GAMS_MODEL.solvestat;

Scalar OBJEST 'best objective', OBJVAL 'objective value';
OBJEST = GAMS_MODEL.objest;
OBJVAL = GAMS_MODEL.objval;

Scalar NUMVAR 'number of variables';
NUMVAR = GAMS_MODEL.numvar

Scalar NUMEQU 'number of equations';
NUMEQU = GAMS_MODEL.numequ

Scalar NUMDVAR 'number of discrete variables';
NUMDVAR = GAMS_MODEL.numdvar

Scalar NUMNZ 'number of nonzeros';
NUMNZ = GAMS_MODEL.numnz

Scalar ETSOLVE 'time to execute solve statement';
ETSOLVE = GAMS_MODEL.etsolve

