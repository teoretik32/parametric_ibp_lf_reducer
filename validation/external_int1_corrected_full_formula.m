(* External Int1 (corrected). The Gamma/EulerGamma prefactor is applied OUTSIDE the
   reducer: PureReduction contains ONLY the certified pure-family reduction of
   Int[(1+x2)^ep*(1+x6)^ep*(1+x2+x6)^(-1+ep), {x2, 0, Infinity}, {x6, 0, Infinity}]. *)
ExternalPrefactor1 = Exp[2*ep*EulerGamma]*Gamma[1-ep]*Gamma[-ep]^2*Gamma[ep]*Gamma[2*ep]/(s*t^2*Gamma[-1-3*ep]*Gamma[1+ep]);

PureReduction = ((4*ep - 1)/(3*(3*ep + 1)))*Int[(x2 + 1)^ep*(x6 + 1)^(ep - 1)*(x2 + x6 + 1)^(ep - 2), {x2, 0, Infinity}, {x6, 0, Infinity}] +
    ((ep - 2)*(5*ep - 2)/(3*ep*(3*ep + 1)))*Int[(x2 + 1)^ep*(x6 + 1)^ep*(x2 + x6 + 1)^(ep - 3), {x2, 0, Infinity}, {x6, 0, Infinity}];

FullIntegralReduction = ExternalPrefactor1*PureReduction;

ReferenceLaurentSeries = 1/(s*t^2)*(1/ep^4 - Pi^2/(12*ep^2) - 43*Zeta[3]/(6*ep) - Pi^4/180);
(* ReferenceLaurentSeries is the expansion around ep=0 (+ O[ep]); reference text
   only -- NOT compared numerically at finite ep. *)
