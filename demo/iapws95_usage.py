"""
 IAPWS-95 Python bindings
"""
from iapws95 import pt2h, pt2s

p, t = 16.10, 535.10

h = pt2h(p, t)
s = pt2s(p, t)
print("(p,t),h,s:",
      "{:>.2f}\t {:>.2f}\t {:>.2f}\t {:>.3f}".format(p, t, h, s))
