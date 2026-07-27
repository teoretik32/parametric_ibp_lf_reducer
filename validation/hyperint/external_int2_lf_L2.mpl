# External Int2 Method.12A -- exact HyperInt input for LF master L2
# label [-1, 0, 0, 0, -1, 0, -1]; certified basis unchanged; audit artifact only.
# REVIEW BEFORE RUNNING: prepared input; do NOT start a long integration
# unattended.  Uncomment the hyperInt call after review.
read "HyperInt.mpl":
_hyper_verbosity := 0:
epsOrder := 3:  # expand through ep^3
f := x2^(ep) * (x2 + 1)^(ep) * (x5 + 1)^(ep - 1) * (x7 + 1)^(-ep - 1) * (r*x2*x5 + x2*x7 + x7 + 1)^(ep - 2):
fser := convert(series(f, ep=0, epsOrder+1), polynom):
intOrder := [x2, x5, x7]:  # valid order from linear-reducibility audit
# result := hyperInt(fser, intOrder):  # <-- uncomment after review
# fibrationBasis(result):
