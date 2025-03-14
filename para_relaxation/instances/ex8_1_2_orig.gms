$offlisting
*  
*  Equation counts
*      Total        E        G        L        N        X        C        B
*          1        1        0        0        0        0        0        0
*  
*  Variable counts
*                   x        b        i      s1s      s2s       sc       si
*      Total     cont   binary  integer     sos1     sos2    scont     sint
*          2        2        0        0        0        0        0        0
*  FX      0
*  
*  Nonzero counts
*      Total    const       NL      DLL
*          2        1        1        0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,objvar;

Positive Variables  x1;

Equations  e1;


e1.. -(588600/POWER(10.8095222429746 - 4.21478541710781*cos((-2.09439333333333)
      + x1),6) - 1079.1/POWER(10.8095222429746 - 4.21478541710781*cos((-
     2.09439333333333) + x1),3) + 600800/POWER(10.8095222429746 - 
     4.21478541710781*cos(x1),6) - 1071.5/POWER(10.8095222429746 - 
     4.21478541710781*cos(x1),3) + 481300/POWER(10.8095222429746 - 
     4.21478541710781*cos(2.09439333333333 + x1),6) - 1064.6/POWER(
     10.8095222429746 - 4.21478541710781*cos(2.09439333333333 + x1),3))
      + objvar =E= 0;

* set non-default bounds
x1.up = 6.28318;

Model m / all /;

m.limrow=0; m.limcol=0;
m.tolproj=0.0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

$if not set NLP $set NLP NLP
Solve m using %NLP% minimizing objvar;
