//
//  euler_spiral.h
//  Euler Spiral
//
//  Created by Sidakpreet Singh Chawla on 13/09/13.
//  Copyright (c) 2013 Sidakpreet Singh Chawla. All rights reserved.
//
#ifndef _EULER_SPIRAL_H_
#define _EULER_SPIRAL_H_

#include "oct/math/math.h"
#include "oct/serv/LogService.h"
#include <cmath>

/* @Brief : This class is used to find the parameters of a euler curve
 * 			between two points given the tangents at the points
 * 			Since there are infinite solutions so the total curvature is
 * 			minimised, this is achieved by keeping the change <= 2 * pi
 * 			Also since no straight away numerical solution exists for the
 * 			euler spiral equations gradient descent is used get the best fit.
 * 			For starting the estimation the initial conditions are to be
 * specified. The following code estimate the initial parameters by fitting a
 * Biarc on the given data. The Curvature in Euler curve is represented as K(s)
 * = gamma*s + k So the purpose is to find k, gamma and the length of the curve
 */

struct EulerSpiral
{
	struct _Private
	{
		int integral_steps;

		// TODO: Replace the following by some calculus algebra libraries for
		// more efficiency
		/* @Brief : Calculates the fresnel integral of x using trapezoidal
		 * 			method of numerical integration
		 * 			Cos term = integration of cos(s^2) with 0 < s < x
		 * 			Sin term = integration of sin(s^2) with 0 < s < x
		 *
		 * @params  [in] x       : The variable whose fresnel integral is to be
		 * found
		 * @params [out] costerm : The value of the cos term of the integration
		 * @params [out] sinterm : The value of the sin term of the integration
		 */
		int fresnel(double x, double & costerm, double & sinterm);

		/* @Brief : Calculates the minimum of the 4 given numbers
		 * 			Returns the minimum of the 4 numbers.
		 * @params  [in] a  : The first number
		 * @params  [in] b  : The Second number
		 * @params  [in] c  : The Third number
		 * @params  [in] d  : The Fourth number
		 */
		double min4(double a, double b, double c, double d);

		/* @Brief : Calculates the minimum of the 3 given numbers
		 * 			Returns the minimum of the 3 numbers.
		 * @params  [in] a  : The first number
		 * @params  [in] b  : The Second number
		 * @params  [in] c  : The Third number
		 */
		double min3(double a, double b, double c);

		// ---------------------------------------------------------------------------
		// 							Biarc Estimator
		// ---------------------------------------------------------------------------

		/* @Brief : This returns the angle b/w two vectors
		 * 			( Sample configuration in complex space - y + ix)
		 * @params  [in] x  : The imaginary part of the vector
		 * @params  [in] y  : The real part of the vector
		 */
		double AngleOfVector(double x, double y);

		/* @Brief : This return the angle at the meeting point of two Biarcs
		 *
		 * @params  [in] x0      : Starting X co-ordinate for biarc formation
		 * @params  [in] y0      : Starting Y co-ordinate for biarc formation
		 * @params  [in] theta0  : The angle of the tangent at the starting
		 * point (Starting orientation)
		 * @params  [in] x2      : Ending X co-ordinate for biarc formation
		 * @params  [in] y2      : Ending Y co-ordinate for biarc formation
		 * @params  [in] theta2  : The angle of the tangent at the Ending point
		 * (End point orientation)
		 * @params  [in] k1      : The curvature at tha starting point
		 * @params  [in] k2      : The curvature at tha ending point
		 */
		double ComputeJoinTheta(double x0,
		                        double y0,
		                        double theta0,
		                        double x2,
		                        double y2,
		                        double theta2,
		                        double k1,
		                        double k2);

		/* @Brief : This returns the arc length for a curve with starting
		 * tagential angle theta0 to the point where the angle of the tangent
		 * vector reaches theta1 with a curve using constant curvature k
		 * @params  [in] theta0  : The angle of the tangent at the starting
		 * point (Starting orientation)
		 * @params  [in] theta1  : The angle of the tangent at the Ending point
		 * (End point orientation)
		 * @params  [in] k       : The curvature of the curve ( Remains
		 * constant)
		 */
		double ComputeArcLength(double theta0, double theta1, double k);

		/* @Brief : Fits a Biarc b/w two points given the tangent vector at both
		 * the points.
		 *
		 * @params  [in] x0        : Starting X co-ordinate for biarc formation
		 * @params  [in] y0        : Starting Y co-ordinate for biarc formation
		 * @params  [in] x2        : Ending X co-ordinate for biarc formation
		 * @params  [in] y2        : Ending Y co-ordinate for biarc formation
		 * @params  [in] theta0    : The angle of the tangent at the starting
		 * point (Starting orientation)
		 * @params  [in] theta2    : The angle of the tangent at the Ending
		 * point (End point orientation)
		 * @params [out] estimateK : The estimate of the curvature from the
		 * Biarc fitting
		 * @params [out] estimatel : The estimate of the length of the Biarc
		 * curve
		 */
		int ComputeBiarcSolution(double x0,
		                         double y0, // Start point coordinates
		                         double x2,
		                         double y2, // End point coordinates
		                         double theta0,
		                         double theta2, // The initial and final angles
		                                        // respectively
		                         double & estimateK,
		                         double & estimatel); // To be calculated

	} _private;
		// ---------------------------------------------------------------------------
		// 						Euler Spiral Estimation
		// ---------------------------------------------------------------------------

		/* @Brief : Function to apply gradient descent for calculation of the
		 * euler curve parameters This takes the initial values from the Biarc
		 * estimation
		 * @params  [in] x0        : Starting X co-ordinate
		 * @params  [in] y0        : Starting Y co-ordinate
		 * @params  [in] theta0    : The angle of the tangent at the starting
		 * point (Starting orientation)
		 * @params  [in] x2        : The required ending X co-ordinate
		 * @params  [in] y2        : The required ending Y co-ordinate
		 * @params  [in] theta2    : The required angle of the tangent at the
		 * Ending point (End point orientation)
		 * @params  [in] estimateK : The estimate of the curvature from the
		 * Biarc fitting
		 * @params  [in] estimatel : The estimate of the length of the Biarc
		 * curve
		 * @params  [in] iternum   : The number of iterations to be done in
		 * gradient descent
		 * @params [out] Kfin      : The curvature parameter of euler curve
		 * @params [out] Lfin      : The length of the euler curve
		 */
		int SolveIteratively(double x0,
		                     double y0,
		                     double theta0,
		                     double x1,
		                     double y1,
		                     double theta2,
		                     double estimateK,
		                     double estimateL,
		                     int iternum,
		                     double & Kfin,
		                     double & Lfin);

		/* @Brief : This function calculates the end points X,Y fromt the curve
		 * paramets k & l and than return the error as the distance as the
		 * straight line distance b/w (X,Y) and (x2,y2)
		 * @params  [in] x0        : Starting X co-ordinate
		 * @params  [in] y0        : Starting Y co-ordinate
		 * @params  [in] theta0    : The angle of the tangent at the starting
		 * point (Starting orientation)
		 * @params  [in] x2        : The required ending X co-ordinate
		 * @params  [in] y2        : The required ending Y co-ordinate
		 * @params  [in] theta2    : The required angle of the tangent at the
		 * Ending point (End point orientation)
		 * @params  [in] k         : The curvature parameter of euler curve
		 * @params  [in] l         : The length of the euler curve
		 */
		double ComputeError(double x0,
		                    double y0,
		                    double theta0,
		                    double x2,
		                    double y2,
		                    double theta2,
		                    double k,
		                    double l);

		/* @Brief : This function calculates the end points X,Y fromt the curve
		 * paramets k & l and than return points X,Y
		 * @params  [in] a         : Starting X co-ordinate
		 * @params  [in] b         : Starting Y co-ordinate
		 * @params  [in] theta     : The angle of the tangent at the starting
		 * point (Starting orientation)
		 * @params  [in] k         : The curvature parameter of euler curve
		 * @params  [in] gamma     : the gamma as in K(s) = gamma*s + k
		 * @params  [in] l         : The length of the euler curve
		 * @params [out] x         : The X co-ordinate of the end point
		 * calculated from the given curve parameters
		 * @params [out] y         : The Y co-ordinate of the end point
		 * calculated from the given curve parameters
		 */
		int EulerSpiralEndPoint(double a,
		                        double b,
		                        double theta,
		                        double k,
		                        double gamma,
		                        double s,
		                        double & x,
		                        double & y);

	EulerSpiral()
	{
		_private.integral_steps = 1000; // Setting the no. of divisions for
		                                // numerical integraion of fresnel
		                                // integration
	}

	/* @Brief : Function to get the euler cuve parameters from the given
	 * starting and end points along with the orientations at both points.
	 * @params  [in] x0        : Starting X co-ordinate
	 * @params  [in] y0        : Starting Y co-ordinate
	 * @params  [in] theta0    : The angle of the tangent at the starting point
	 * (Starting orientation)
	 * @params  [in] x2        : The required ending X co-ordinate
	 * @params  [in] y2        : The required ending Y co-ordinate
	 * @params  [in] theta2    : The required angle of the tangent at the Ending
	 * point (End point orientation)
	 * @params  [in] iternum   : The number of iterations to be done in gradient
	 * descent
	 * @params [out] Kfin      : The curvature parameter of euler curve
	 * @params [out] Lfin      : The length of the euler curve
	 */
	int SolveEuler(double x0,
	               double y0,
	               double theta0,
	               double x1,
	               double y1,
	               double theta2,
	               int iternum,
	               double & Kfin,
	               double & Lfin,
	               double & Gfin);

	/* @Brief : Function to set the no. of divisions for numerical integraion of
	 * fresnel integration This is required as for large difference b/w (x2,y2)
	 * and (x0,y0) the default value 1000 is not enough
	 * @params  [in] no_of_integral_steps    : No. of divisions for numerical
	 * integraion of fresnel integration
	 */
	void setFresenlIntegrationSteps(int no_of_integral_steps);
};

#endif
