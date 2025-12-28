* HW8: Fatih ŞAHİNOĞLU
* Fit S = a*(h-h0)^b by least squares (relative heights)

SET i / i1*i12 /;

SCALAR h0 "reference level" / 0 /;

PARAMETER
    h(i)        "Height (relative units)"
    Sderived(i) "Storage derived from QGIS (Mm3)";

* Data (decimal comma converted to decimal dot)
h('i1')  = 0;  Sderived('i1')  = 0.000;
h('i2')  = 5;  Sderived('i2')  = 79.600;
h('i3')  = 10;  Sderived('i3')  = 151.300;
h('i4')  = 15;  Sderived('i4')  = 216.400;
h('i5')  = 20;  Sderived('i5')  = 272.400;
h('i6')  = 25;  Sderived('i6')  = 311.400;
h('i7')  = 30;  Sderived('i7')  = 340.300;
h('i8')  = 35;  Sderived('i8')  = 366.900;
h('i9')  = 40;  Sderived('i9')  = 396.400;
h('i10') = 45;  Sderived('i10') = 414.600;
h('i11') = 50;  Sderived('i11') = 436.600;
h('i12') = 55;  Sderived('i12') = 478.700;

POSITIVE VARIABLES
    a "coefficient"
    b "exponent";

VARIABLES
    e(i) "residual"
    obj  "sum of squared residuals";

EQUATIONS
    residual(i)
    objective;

* Only fit points with (h-h0) > 0 => exclude i1 (where h=h0)
residual(i)$(ord(i) > 1)..
    e(i) =E= a * exp( b * log(h(i) - h0) ) - Sderived(i);

objective..
    obj =E= sum(i$(ord(i) > 1), sqr(e(i)));

MODEL fitSE / residual, objective /;

* Starting values (just initial guesses)
a.L = 5;
b.L = 2;

SOLVE fitSE USING NLP MINIMIZING obj;

DISPLAY a.L, b.L, obj.L;

* OPTIONAL: Export results for plotting (CSV)
FILE out / "HW8_fit_output.csv" /;
PUT out;
PUT "h,Sderived,Sfit" /;

LOOP(i,
    IF(ord(i) = 1,
        PUT h(i):0:6 "," Sderived(i):0:6 "," 0:0:6 /;
    ELSE
        PUT h(i):0:6 "," Sderived(i):0:6 "," (a.L * exp(b.L * log(h(i) - h0))):0:6 /;
    );
);

PUTCLOSE out;
