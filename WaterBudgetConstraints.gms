* ------------------------------------------------------------
* Reservoir Operation - Year 2000 (Constrained)
* Release differ from Demand when operations become constrained.

SET t /1*12/;

SCALAR
    K      "Max storage capacity (Mm3)" /478.7/,
    Smin0  "Min storage (Mm3)"          /0/,
    m      "Slope A(S)=m*S (km2/Mm3)"   /0.0363633/,
    beg_s  "Initial storage (Mm3) set as 0.5*K",
    fracHi "Normal max release fraction of available storage" /0.06/,
    fracLo "Low-storage max release fraction of available storage" /0.04/,
    thr    "Threshold = 50% of initial storage (Mm3)",
    ksig   "Sigmoid steepness (Mm3) (smaller => sharper switch)" /5/;

beg_s = 0.5*K;
thr   = 0.5*beg_s;

PARAMETER
    Q(t)       "Monthly inflow (Mm3/month)",
    D(t)       "Monthly demand (Mm3/month)",
    Evap_mm(t) "Monthly evaporation depth (mm/month)",
    a(t)       "Evap coefficient",
    b(t)       "Evap constant term (Mm3/month)";

* Inflow Q (Mm3/month)
Q('1')  = 3.292903;
Q('2')  = 4.880690;
Q('3')  = 9.136129;
Q('4')  = 14.828667;
Q('5')  = 11.873226;
Q('6')  = 8.896118;
Q('7')  = 5.905753;
Q('8')  = 4.365127;
Q('9')  = 3.688102;
Q('10') = 3.333168;
Q('11') = 2.981884;
Q('12') = 2.879407;

* Demand (constant)
D(t) = 14;

* Evaporation depth (mm/month)
Evap_mm('1')  = 30;
Evap_mm('2')  = 35;
Evap_mm('3')  = 80;
Evap_mm('4')  = 100;
Evap_mm('5')  = 112;
Evap_mm('6')  = 160;
Evap_mm('7')  = 190;
Evap_mm('8')  = 175;
Evap_mm('9')  = 90;
Evap_mm('10') = 80;
Evap_mm('11') = 63;
Evap_mm('12') = 51;

a(t) = 0.5 * (Evap_mm(t)/1000) * m;
b(t) = 0;

POSITIVE VARIABLES
    S(t)      "End of month storage (Mm3)",
    R(t)      "Release supplied (Mm3/month)",
    Short(t)  "Shortage (Mm3/month)",
    Savail(t) "Available storage at start of month (Mm3)";

VARIABLE
    z         "Objective";

EQUATIONS
    obj,
    balance(t),
    availDef(t),
    demandBal(t),
    relCap(t);

obj.. z =e= sum(t, sqr(Short(t)));

availDef(t)..
    Savail(t) =e= beg_s$(ord(t)=1) + S(t-1)$(ord(t)>1);

demandBal(t)..  R(t) + Short(t) =e= D(t);

relCap(t)..
    R(t) =l=
      ( fracLo
        + (fracHi - fracLo) * ( 1 / (1 + exp(-(Savail(t) - thr)/ksig)) )
      ) * Savail(t);

balance(t)..
    (1 + a(t)) * S(t) =e=
        (1 - a(t)) * Savail(t)
        + Q(t) - R(t) - b(t);

* Bounds
S.lo(t) = Smin0;
S.up(t) = K;
R.lo(t) = 0;

MODEL ResOpPolicy /all/;
SOLVE ResOpPolicy USING NLP MINIMIZING z;

* Storage change for plotting
PARAMETER dS(t) "Change in storage (Mm3/month)";
dS(t) = S.l(t) - Savail.l(t);

DISPLAY K, beg_s, thr, fracHi, fracLo, ksig, Q, D, Evap_mm, a,
        Savail.l, R.l, Short.l, S.l, dS, z.l;

* Clean CSV export
FILE fout /resop_results_2000_constrained_policy.csv/;
PUT fout;
PUT "Month,Q,D,Evap_mm,a,Savail,R,Short,S_end,dS" /;

LOOP(t,
  PUT ord(t), ",",
      Q(t), ",",
      D(t), ",",
      Evap_mm(t), ",",
      a(t), ",",
      Savail.l(t), ",",
      R.l(t), ",",
      Short.l(t), ",",
      S.l(t), ",",
      dS(t) /
);

PUTCLOSE fout;
