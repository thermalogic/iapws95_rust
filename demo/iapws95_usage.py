"""
 IAPWS-95 Python bindings
"""
from iapws95 import  tr2h, tr2p

t,r = 300-273.15,0.9965560e3

p = tr2p(t,r)
h = tr2h(t,r)

print("(t,r),p,h:",
      "{:>.2f}\t {:>.2f}\t {:>.2f}\t {:>.3f}".format(t,r, p,h))
