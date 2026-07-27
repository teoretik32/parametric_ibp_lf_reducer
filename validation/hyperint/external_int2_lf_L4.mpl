# External Int2 Method.12A -- exact HyperInt input for LF master L4
# label [0, 0, 1, -1, 0, 0, -1]; certified basis unchanged; audit artifact only.
# REVIEW BEFORE RUNNING: prepared input; do NOT start a long integration
# unattended.  Uncomment the hyperInt call after review.
read "HyperInt.mpl":
_hyper_verbosity := 0:
epsOrder := 4:  # expand through ep^4
f := x2^(ep + 1) * x7^(1) * (x2 + 1)^(ep - 1) * (x5 + 1)^(ep) * (x7 + 1)^(-ep - 1) * (r*x2*x5 + x2*x7 + x7 + 1)^(ep - 2):
fser := convert(series(f, ep=0, epsOrder+1), polynom):
intOrder := [x2, x5, x7]:  # valid order from linear-reducibility audit
# result := hyperInt(fser, intOrder):  # <-- uncomment after review
# fibrationBasis(result):
