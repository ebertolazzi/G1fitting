### G1 fitting with clothoids
**by Enrico Bertolazzi and Marco Frego**

The script buildClothoid implements the algorithm described in the paper

*G1 fitting with clothoids*, Mathematical Methods in the Applied Sciences, John Wiley & Sons, (2014), Ltd,.
http://onlinelibrary.wiley.com/doi/10.1002/mma.3114/abstract


**Description:**
Given two points and two direction associated with the points, 
a clothoid, i.e. a curve with linear varying curvature is computed
in such a way it pass to the points with the prescribed direction.
The solution in general is not unique but chosing the one for
which the angle direction variation is less than `2*pi` the solution
is unique.
The sofware solves the nonlinear system associated to the fitting problem
computing initial curvature and its derivative with the lenght of the curve.
An additional routine for the computation of the points along a clothoid
curve is added for convenience.

**Usage:**
To compute curvature and length of a clothoid passing throught points
P0=`(x0,y0)`, P1=`(x1,y1)` with angles `theta0` and `theta1` use the
function `buildClothoid`

`[k,dk,L] = buildClothoid( x0, y0, theta0, x1, y1, theta1, tol ) ;`

The parameter `tol` (usually `1e-10`) is a tolerance parameter
used to stop Newton iteration.
The resulting curve can be described by the 5 parameters

  - `(x0,y0)` initial point
  - `theta0`  initial direction (angle)
  - `k`       initial curvature of the curve
  - `dk`      derivative of the cuvature along arc length

plus a 6th parameter `L`

  - L         total length of the curve connecting P0 and P1

to compute points along a clothoid curve use the function `pointsOnClothoid`

`XY = pointsOnClothoid( x0, y0, theta0, k, dk, L, npts ) ;`

This function uses the 5 parameters `x0`, `y0`, `theta0`, `k`, `dk`
which indentify the curve. The parameter `L` is used to determine length 
of the portion of the curve to compute. The parameter `npts` is the number
of points computed along the curve. 

`XY` is a `2 x npts` matrix whose columns are the points along the curve. 

To plot the computed curve use MATLAB `plot` command as usual:

`plot( XY(1,:), XY(2,:), '-r' ) ;`

Three sample scripts: TestN0, TestN1, TestN2 shows how to use the functions.

**Authors:**
	
	Enrico Bertolazzi and Marco Frego
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
	m.fregox@gmail.com

