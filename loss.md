# TD2 Loss Coefficients

This page demonstrates how TD2 computes loss. TD2 uses two variables to hold the loss values `YCON` and `YOSS` according to the original manual.

**YCON** (Array[Float]): Value of the 9 constants which define the internal loss correlation. Example: [0.04, 0.0, 0.5, 1.0, 1.5, 0.04, 0.0, 1.0, 2]

if $V_in/V_ex$ \geq a_3$
$$ Y = \frac{|tan(\beta_{in}) - tan(\beta_{ex})|}{a_4+a_5 cos(\beta_ex)} \bigg[a_1 + a_2(V_{in}/V_{ex}-a_3)\bigg]\frac{w}{2E6 \mu r}^{-0.2}
$$

else $V_in/V_ex$ < a_3$

$$ Y = \frac{|tan(\beta_{in}) - tan(\beta_{ex})|}{a_4 + a_5 cos(\beta_{ex})} \bigg[ a_6 + a_7 \frac{V_{in}}{V_{ex}}^a_8 \bigg] \frac{w}{2E6 \mu r}^{-0.2}
$$

Recommended values are:
a_1 = a_6 = 0.057
a_2 = a_7 = 0
a_3 = a_8 = 0 or arbitrary
a_4 = 1.0
a_5 = 1.5

**YOSS** (Array[Array[Float]]): Values of the loss coefficient (if ISPEC=0) or an **additional loss factor** (if ISPEC=2) at each stator exit of the spool corresponding to the radial
coordinates RNXT. This is a modifier for the pressure loss. In the original code YOSS -> YOS -> FACL which is then used to modify the Pressure loss.

td2-2.f line 557
`YOS(J)=YOSS(J,IBR)`


td2-2.f line 2765
`WYECOR = FACL(I)*WYECOR(I)`


td2-2.f line 1415

```fortran
CALL I1AP1(RP,YOSP,RXTS,YOS,NXT) ! Interpolation of YOS
FACL(J)=YOSP ! YOSP is the interpolated values 
```

