"""
 IAPWS-95 Python bindings
"""
from iapws95 import pt2h, pt2s, tr2h, tr2s, tr2cv, tr2cp, tr2w, tr2x

p, t = 16.10, 535.10

h = pt2h(p, t)
s = pt2s(p, t)
print("(p,t),h,s:",
      "{:>.2f}\t {:>.2f}\t {:>.2f}\t {:>.3f}".format(p, t, h, s))
