# edgecompare

This is a search program for edge lengths of hybrid Archimedean tilings of the hyperbolic plane.  The algorithm, in which vertex configurations are checked in ascending order of edge length, was conceived by Marek, who is active on the HyperRogue Discord.  I have translated it into Rust for better performance.

## Background
On the Euclidean plane, we can make tessellations by putting 6 regular triangles around each vertex, or by putting 4 squares around each vertex, or 3 regular hexagons.  If we put even more regular polygons around a vertex, this may form a [tiling of the hyperbolic plane](https://en.wikipedia.org/wiki/Template:Regular_hyperbolic_tiling_table).  The hyperbolic plane has constant negative curvature. For a specific configuration of regular polygons to tile the hyperbolic plane, the polygons' edge length must take on a specific value.

For a tiling with $q$ regular $p$-gons around each vertex, the side length $s$ of the polygons must be
$$ s = 2 \mathrm{arcosh}\left( \frac{\cos \frac{\pi}{p} }{\sin \frac{\pi}{q} } \right). $$

The interior angles $\theta = 2\pi/q$ must satisfy
$$\theta = 2 \arcsin \left( \frac{\cos \frac{\pi}{p}}{\cosh \frac{s}{2}} \right).$$

Now consider some other tiling of the hyperbolic plane, where we fit $N_g$ groups of regular polygons around a vertex, and the $i$-th group consists of $q_i$ regular $p_i$-gons with side length $k_is$, where $k_i$ is a whole number.  Then, since the angles around a circle must still add up to $2\pi$ radians, the base side length $s$ is determined by
$$\sum_{i=1}^{N_g} 2 q_i \arcsin \left( \frac{\cos \frac{\pi}{p_i}}{\cosh \frac{k_is}{2}} \right) = 2\pi.$$

We allow for $p_i$ to be infinite - in this case, $\cos \frac{\pi}{p_i} \rightarrow 1$.  This is called an apeirogon.

In this program, given a vertex configuration $(k_i, p_i, q_i)_{1 \leq i \leq N_g}$, we solve for $s$ numerically.  We enumerate many vertex configurations and calculate $s$ for each of them, looking for values of $s$ that correspond to multiple vertex configurations, since those could lead to hybrid tilings.

For convenience, we define $m = \cosh(s/2)$.