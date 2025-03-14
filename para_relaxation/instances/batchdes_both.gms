$offlisting
$offdigit

EQUATIONS
	c1_lo
	c2_lo
	c3_lo
	c4_lo
	c5_lo
	c6_lo
	c7_lo
	c8_lo
	c9_lo
	c10_lo
	c11_lo
	c12_lo
	c13_hi
	c14
	c15
	c16
	c17
	c18
	c19
	c20
	c21
	c22
	c23
	c24
	c25
	c26
	c27
	c28
	c29
	c30
	c31
	c32
	c33
	c34
	c35_hi
	c36_hi
	c37_lo
	c38_lo
	c39_hi
	c40_hi
	c41_lo
	c42_lo
	c43;

BINARY VARIABLES
	x13
	x14
	x15
	x16
	x17
	x18
	x19
	x20
	x21;

POSITIVE VARIABLES
	x7
	x8
	x9;

VARIABLES
	GAMS_OBJECTIVE
	x1
	x2
	x3
	x4
	x5
	x6
	x10
	x11
	x12
	x22
	x23
	x24
	x25
	x26
	x27
	x28
	x29
	x30
	x31
	x32
	x33
	x34;


c1_lo.. 0.6931471805599449 =l= x1 - x2 ;
c2_lo.. 1.09861228866811 =l= x3 - x2 ;
c3_lo.. 1.3862943611198899 =l= x4 - x2 ;
c4_lo.. 1.3862943611198899 =l= x1 - x5 ;
c5_lo.. 1.79175946922805 =l= x3 - x5 ;
c6_lo.. 1.09861228866811 =l= x4 - x5 ;
c7_lo.. 2.07944154167984 =l= x6 + x7 ;
c8_lo.. 2.99573227355399 =l= x6 + x8 ;
c9_lo.. 1.3862943611198899 =l= x6 + x9 ;
c10_lo.. 2.3025850929940499 =l= x10 + x7 ;
c11_lo.. 2.4849066497879999 =l= x10 + x8 ;
c12_lo.. 1.09861228866811 =l= x10 + x9 ;
c13_hi.. x11 + x12 =l= 6000 ;
c14.. (-0.6931471805599449)*x13 + (-1.09861228866811)*x14 + x7 =e= 0 ;
c15.. (-0.6931471805599449)*x15 + (-1.09861228866811)*x16 + x8 =e= 0 ;
c16.. (-0.6931471805599449)*x17 + (-1.09861228866811)*x18 + x9 =e= 0 ;
c17.. x19 + x13 + x14 =e= 1 ;
c18.. x20 + x15 + x16 =e= 1 ;
c19.. x21 + x17 + x18 =e= 1 ;
c20.. -x22 + 250*x23 =e= 0 ;
c21.. -x24 + 500*x25 =e= 0 ;
c22.. -x26 + 340*x27 =e= 0 ;
c23.. -x11 + 200000*x28 =e= 0 ;
c24.. -x12 + 150000*x29 =e= 0 ;
c25.. -x23 + exp(x30) =e= 0 ;
c26.. -x25 + exp(x31) =e= 0 ;
c27.. -x27 + exp(x32) =e= 0 ;
c28.. -x28 + exp(x33) =e= 0 ;
c29.. -x29 + exp(x34) =e= 0 ;
c30.. 0.59999999999999998*x1 + x7 - x30 =e= 0 ;
c31.. 0.59999999999999998*x3 + x8 - x31 =e= 0 ;
c32.. 0.59999999999999998*x4 + x9 - x32 =e= 0 ;
c33.. -x2 + x6 - x33 =e= 0 ;
c34.. -x5 + x10 - x34 =e= 0 ;
c35_hi.. 0.15093999999999999*x33 - x28 + 0.016199999999999999*x33*x33 =l= -0.35642694470259678 ;
c36_hi.. 0.18733*x33 - x28 + 0.020959999999999999*x33*x33 =l= -0.41938694669644855 ;
c37_lo.. -0.39146428517606907 =l= 0.163*x33 - x28 + 0.017469999999999999*x33*x33 ;
c38_lo.. -0.66975564820580558 =l= 0.37211*x33 - x28 + 0.05523*x33*x33 ;
c39_hi.. 0.15093999999999999*x34 - x29 + 0.016199999999999999*x34*x34 =l= -0.35642694470259678 ;
c40_hi.. 0.18733*x34 - x29 + 0.020959999999999999*x34*x34 =l= -0.41938694669644855 ;
c41_lo.. -0.39146428517606907 =l= 0.163*x34 - x29 + 0.017469999999999999*x34*x34 ;
c42_lo.. -0.66975564820580558 =l= 0.37211*x34 - x29 + 0.05523*x34*x34 ;
c43.. GAMS_OBJECTIVE =e= x22 + x24 + x26 ;

x7.up = 1.09861228866811;
x8.up = 1.09861228866811;
x9.up = 1.09861228866811;
x1.lo = 5.5214609178622496;
x1.up = 7.82404601085629;
x2.lo = 5.40367788220586;
x2.up = 6.4377516497364;
x3.lo = 5.5214609178622496;
x3.up = 7.82404601085629;
x4.lo = 5.5214609178622496;
x4.up = 7.82404601085629;
x5.lo = 4.60517018598809;
x5.up = 6.03228654162824;
x6.lo = 1.89711998488588;
x6.up = 2.99573227355399;
x10.lo = 1.3862943611198899;
x10.up = 2.4849066497879999;
x11.lo = 2133.333333333332;
x11.up = 18000.000000000025;
x12.lo = 1439.9999999999948;
x12.up = 18000.000000000015;
x22.lo = 6866.0033956632478;
x22.up = 82002.15554574574;
x23.lo = 27.464013582652992;
x23.up = 328.00862218298295;
x24.lo = 13732.006791326496;
x24.up = 164004.31109149149;
x25.lo = 27.464013582652992;
x25.up = 328.00862218298295;
x26.lo = 9337.764618102017;
x26.up = 111522.9315422142;
x27.lo = 27.464013582652992;
x27.up = 328.00862218298295;
x28.lo = 0.010666666666666661;
x28.up = 0.09000000000000013;
x29.lo = 0.009599999999999964;
x29.up = 0.1200000000000001;
x30.lo = 3.3128765507173497;
x30.up = 5.7930398951818836;
x31.lo = 3.3128765507173497;
x31.up = 5.7930398951818836;
x32.lo = 3.3128765507173497;
x32.up = 5.7930398951818836;
x33.lo = -4.5406316648505207;
x33.up = -2.4079456086518705;
x34.lo = -4.64599218050835;
x34.up = -2.12026353620009;

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

