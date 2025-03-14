$offlisting
*  
*  Equation counts
*      Total        E        G        L        N        X        C        B
*          9        9        0        0        0        0        0        0
*  
*  Variable counts
*                   x        b        i      s1s      s2s       sc       si
*      Total     cont   binary  integer     sos1     sos2    scont     sint
*         15       15        0        0        0        0        0        0
*  FX      0
*  
*  Nonzero counts
*      Total    const       NL      DLL
*         65        9       56        0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,objvar;

Positive Variables  x1,x2,x3,x4,x5,x6,x7,x8,x12,x13,x14;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9;


e1.. -(sqr((-0.1622 + x1)/x1) + sqr((-0.6791 + x2)/x2) + sqr((-0.679 + x3)/x3)
      + sqr((-0.3875 + x4)/x4) + sqr((-0.1822 + x5)/x5) + sqr((-0.1249 + x6)/x6
     ) + sqr((-0.0857 + x7)/x7) + sqr((-0.0616 + x8)/x8)) + objvar =E= 0;

e2.. exp(-4*x12)*x9 + exp(-4*x13)*x10 + exp(-4*x14)*x11 - x1 =E= 0;

e3.. exp(-8*x12)*x9 + exp(-8*x13)*x10 + exp(-8*x14)*x11 - x2 =E= 0;

e4.. exp(-12*x12)*x9 + exp(-12*x13)*x10 + exp(-12*x14)*x11 - x3 =E= 0;

e5.. exp(-24*x12)*x9 + exp(-24*x13)*x10 + exp(-24*x14)*x11 - x4 =E= 0;

e6.. exp(-48*x12)*x9 + exp(-48*x13)*x10 + exp(-48*x14)*x11 - x5 =E= 0;

e7.. exp(-72*x12)*x9 + exp(-72*x13)*x10 + exp(-72*x14)*x11 - x6 =E= 0;

e8.. exp(-94*x12)*x9 + exp(-94*x13)*x10 + exp(-94*x14)*x11 - x7 =E= 0;

e9.. exp(-118*x12)*x9 + exp(-118*x13)*x10 + exp(-118*x14)*x11 - x8 =E= 0;

* set non-default bounds
x1.up = 1;
x2.up = 1;
x3.up = 1;
x4.up = 1;
x5.up = 1;
x6.up = 1;
x7.up = 1;
x8.up = 1;
x9.lo = -10; x9.up = 10;
x10.lo = -10; x10.up = 10;
x11.lo = -10; x11.up = 10;
x12.up = 0.5;
x13.up = 0.5;
x14.up = 0.5;

* set non-default levels
x1.l = 0.171747132;
x2.l = 0.843266708;
x3.l = 0.550375356;
x4.l = 0.301137904;
x5.l = 0.292212117;
x6.l = 0.224052867;
x7.l = 0.349830504;
x8.l = 0.856270347;
x9.l = 0.355;
x10.l = 2.007;
x11.l = -4.575;
x12.l = 0.015;
x13.l = 0.11;
x14.l = 0.285;

Model m / all /;

m.limrow=0; m.limcol=0;
m.tolproj=0.0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

$if not set NLP $set NLP NLP
Solve m using %NLP% minimizing objvar;
