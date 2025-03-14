$offlisting
*  
*  Equation counts
*      Total        E        G        L        N        X        C        B
*          1        1        0        0        0        0        0        0
*  
*  Variable counts
*                   x        b        i      s1s      s2s       sc       si
*      Total     cont   binary  integer     sos1     sos2    scont     sint
*          3        3        0        0        0        0        0        0
*  FX      0
*  
*  Nonzero counts
*      Total    const       NL      DLL
*          3        1        2        0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,objvar;

Equations  e1;


e1.. -(exp(sin(50*x1)) + sin(60*exp(x2)) + sin(70*sin(x1)) + sin(sin(80*x2)) - 
     sin(10*x1 + 10*x2) + 0.25*(sqr(x1) + sqr(x2))) + objvar =E= 0;

* set non-default bounds
x1.lo = -3; x1.up = 3;
x2.lo = -3; x2.up = 3;

* set non-default levels
x1.l = -0.655668942;
x2.l = 0.346914252;

Model m / all /;

m.limrow=0; m.limcol=0;
m.tolproj=0.0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

$if not set NLP $set NLP NLP
Solve m using %NLP% minimizing objvar;
