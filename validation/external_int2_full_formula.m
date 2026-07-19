(* External Int2 (dimensionless, r = s/t). The prefactor P2 (Gamma ratio and the
   t^(-3-ep) scaling) is applied OUTSIDE the reducer: PureReduction contains ONLY
   the certified pure-family reduction of Int[F2, x2>0, x5>0, x7>0] with
   F2 = x2^(1+ep)*(1+x2)^ep*(1+x5)^ep*(1+x7)^(-1-ep)*(1+x7+x2*x7+r*x2*x5)^(-1+ep). *)
ExternalPrefactor2 = t^(-3-ep)*Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/(Gamma[-1-3*ep]*Gamma[-2*ep]);

PureReduction = ((ep - 1)/(6*ep*r))*Int[x2^ep*(x2 + 1)^(ep - 1)*(x5 + 1)^(ep - 1)*(x7 + 1)^(-ep - 1)*(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 1), {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}] +
    (-(ep - 1)^2*(r + 1)/(6*ep^2*r))*Int[x2^ep*(x2 + 1)^ep*(x5 + 1)^(ep - 1)*(x7 + 1)^(-ep - 1)*(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 2), {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}] +
    ((ep - 1)*(2*ep^2*r + ep^2 - 1)/(6*ep^3*r))*Int[x2^(ep + 1)*(x2 + 1)^(ep - 1)*(x5 + 1)^ep*(x7 + 1)^(-ep - 1)*(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 2), {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}] +
    (1)*Int[x2^(ep + 1)*(x2 + 1)^ep*(x5 + 1)^(ep - 1)*(x7 + 1)^(-ep - 1)*(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 1), {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}] +
    (-(ep + 1)/ep)*Int[x2^(ep + 1)*(x2 + 1)^ep*(x5 + 1)^ep*(x7 + 1)^(-ep - 2)*(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 1), {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}] +
    ((ep - 1)*(ep + 1)*(2*ep - 1)/(6*ep^3*r))*Int[x2^(ep + 1)*x7*(x2 + 1)^(ep - 1)*(x5 + 1)^ep*(x7 + 1)^(-ep - 1)*(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 2), {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}];

FullIntegralReduction = ExternalPrefactor2*PureReduction;

(* Reference value: AnsvInt2 is not available in this repository and is NOT
   invented here. If a source reference (GPL G[...] expression) is added later it
   must stay metadata only -- GPL values are never reducer coefficients. *)
