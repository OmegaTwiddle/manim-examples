import random
import numpy as np
from math import e
from manim import *

config.frame_width=40

PLANE_FACTOR=25
PLANE_AXES=[-4,4]
COLOR_ORDER=[RED,GRAY,ORANGE,BLUE,GREEN]

POLY_SIZE=17

def angle(c):
    return (np.angle(c)) % (2*PI)

def colors(i):
    return COLOR_ORDER[i%len(COLOR_ORDER)]
def grays(i):
    return WHITE

def get_nth_roots(n):
    nth_poly = [1] + [0]*(n-1) + [-1]
    return sorted(np.roots(nth_poly), key=angle)

def stitch(p1, p2):
    res = []
    for i in range(0,len(p1)):
        res.append(p1[i])
        res.append(p2[i])
        res.append(p1[(i+1)%len(p1)])
    return res
        

def gen_new_poly(points, nthroot, kthpow):
    res = []
    nthroots = get_nth_roots(nthroot)
    a = nthroots[1] ** kthpow
    
    for i in range(0,len(points)):
        p1 = points[i]
        p2 = points[(i+1)%len(points)]
        new_p = ((1-a)**(-1)) * (p2 - a*p1)
        res.append(new_p)
    return res
    

def gen_random_curve(n):
    theta = 0
    points = []
    for i in range(0,n):
        theta += (2 * PI) / n 
        mag = 1 + 0.2 * (0.5-random.random())
        p = mag * e**((0+1j)*theta)
        points.append(p)
    return points


class PetrExample2(Scene):
    def build_vgroup(self, plane, points, colors_fn, op=1):
        poly = VGroup()
        for i in range(0,len(points)):
            r1 = points[i]
            r2 = points[(i+1)%len(points)]
            l = Line(plane.n2p(r1),plane.n2p(r2))
            l.set_opacity(op)
            d = Dot(plane.n2p(r1)).set_color(colors_fn(i))
            d.set_opacity(op)
            poly.add(l,d)
        return [poly,points]

    def construct(self):
        self.clear()
        plane = ComplexPlane(
            x_range=PLANE_AXES,
            y_range=PLANE_AXES,
            x_length=PLANE_FACTOR,
            y_length=PLANE_FACTOR,
        ).add_coordinates()
        self.add(plane)

        #poly=[0,1,1-1j,0-1j,-1+1j]
        poly = gen_random_curve(POLY_SIZE)
        [vg1,p1] = self.build_vgroup(plane, poly, colors)
        self.add(vg1)
        self.wait(2)

        vg_old = vg1
        orig = range(1,POLY_SIZE-2 + 1)
        for i in reversed(orig):
            # Create new polygon
            new_poly = gen_new_poly(poly, POLY_SIZE, i)
            stitched_poly = stitch(poly,new_poly)

            # Create stitched polygon
            [stitched_vg,p2] = self.build_vgroup(plane, stitched_poly, grays, op=0.3)

            # Create next polygon
            [vg3,p3] = self.build_vgroup(plane, new_poly, colors)

            # Animate
            RT=0.5
            RT2=0.4
            self.play(FadeIn(stitched_vg, run_time=RT))
            self.play(FadeIn(vg3))
            self.play(
                FadeOut(vg_old, run_time=RT),
                FadeOut(stitched_vg, run_time=RT2),
            )
            
            # rescale if needed
            centroid = np.mean(new_poly)
            scalefac = max([np.abs(x-centroid) for x in new_poly])
            if scalefac < 0.7 or scalefac > 7:
                scaled_poly = [(1/scalefac) * (x-centroid) for x in new_poly]
                [vg3_new,pnew] = self.build_vgroup(plane, scaled_poly, colors)
                self.play(
                    ReplacementTransform(vg3, vg3_new, lag_ratio=0, run_time=0.2, rate_func=linear),
                )
                new_poly = pnew
                vg3 = vg3_new

            # Prepare for tomorrow
            poly = new_poly
            vg_old = vg3
        self.wait(4)

