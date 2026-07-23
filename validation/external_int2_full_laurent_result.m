(* External Int2 compact Laurent result through ep^0.
   r = s/t; H[word] denotes G[word; r] with alphabet {0,-1}. *)

r = s/t;
L = Log[r];
T = Log[t];

Fminus4 = -4;
Fminus3 = 3*L + 4*T;
Fminus2 = 5*Pi^2/3 - L^2/2 - 3*L*T - 2*T^2;

Fminus1 =
    2*G[-1,0,0,r]
  - L^3/2 + L^2*T/2 + 3*L*T^2/2 + 2*T^3/3
  - 7*Pi^2*L/3 - 5*Pi^2*T/3
  + Pi^2*G[-1,r] + 62*Zeta[3]/3;

F0 =
    11*L^4/24 + L^3*T/2
  + L^2*(-T^2/4 + 11*Pi^2/6)
  + L*(-T^3/2 + 7*Pi^2*T/3 - 18*Zeta[3])
  - T^4/6 + 5*Pi^2*T^2/6
  - 2*T*G[-1,0,0,r] - 62*T*Zeta[3]/3
  - 4*Pi^2*G[0,-1,r] - 8*G[0,-1,0,0,r]
  - 4*Pi^2*G[-1,0,r]/3 - 6*G[-1,0,0,0,r]
  + Pi^2*G[-1,-1,r] + 2*G[-1,-1,0,0,r]
  + G[-1,r]*(-Pi^2*T + 2*Zeta[3])
  + 23*Pi^4/45;

ExternalInt2Laurent =
  1/(s*t^2) * (
      Fminus4/ep^4
    + Fminus3/ep^3
    + Fminus2/ep^2
    + Fminus1/ep
    + F0
    + O[ep]
  );
