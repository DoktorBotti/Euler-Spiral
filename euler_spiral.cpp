#include "euler_spiral.h"
int EulerSpiral::_Private::fresnel(double x, double & costerm, double & sinterm)
{
	costerm = 0;
	sinterm = 0;
	costerm += (cos(0) + cos((x * x))) / 2;
	sinterm += (sin(0) + sin((x * x))) / 2;

	for (int i = 0; i < integral_steps; i++)
	{
		costerm += cos(((x)*i / integral_steps) * ((x)*i / integral_steps));
		sinterm += sin(((x)*i / integral_steps) * ((x)*i / integral_steps));
	}

	costerm = costerm * (x) / integral_steps;
	sinterm = sinterm * (x) / integral_steps;
	return 0;
}

double EulerSpiral::_Private::min4(double a, double b, double c, double d)
{
	if (a <= b && a <= c && a <= d)
		return a;

	if (b <= a && b <= c && b <= d)
		return b;

	if (c <= a && c <= b && c <= d)
		return c;

	return d;
}

double EulerSpiral::_Private::min3(double a, double b, double c)
{
	if (a <= b && a <= c)
		return a;

	if (b <= a && b <= c)
		return b;

	return c;
}

double EulerSpiral::_Private::AngleOfVector(double x, double y)
{
	double angle = atan2(y, x);
	if (angle < 0)
	{
		angle = angle + 2 * oct::math::pi;
	}
	return angle;
}

double EulerSpiral::_Private::ComputeJoinTheta(double x0,
                                               double y0,
                                               double theta0,
                                               double x2,
                                               double y2,
                                               double theta2,
                                               double k1,
                                               double k2)
{
	double sin_thetaJoin =
	    (k1 * k2 * (x2 - x0) + k2 * sin(theta0) - k1 * sin(theta2)) / (k2 - k1);
	double cos_thetaJoin =
	    (k1 * k2 * (y2 - y0) + k2 * cos(theta0) - k1 * cos(theta2)) / (k2 - k1);

	return AngleOfVector(sin_thetaJoin, cos_thetaJoin);
}

double
EulerSpiral::_Private::ComputeArcLength(double theta0, double theta1, double k)
{
	double numerator = theta1 - theta0;
	if (k < 0 && numerator > 0)
	{
		numerator = numerator - 2 * oct::math::pi;
	}
	else if (k > 0 && numerator < 0)
	{
		numerator = numerator + 2 * oct::math::pi;
	}
	return numerator / k;
}

int EulerSpiral::_Private::ComputeBiarcSolution(
    double x0,
    double y0, // Start point coordinates
    double x2,
    double y2, // End point coordinates
    double theta0,
    double theta2, // The initial and final angles respectively
    double & estimateK,
    double & estimatel) // To be calculated

{
	double L = sqrt((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0));
	double psi = AngleOfVector(x2 - x0, y2 - y0);

	double smallAngle = 0.001;
	double smallCurvature = 0.001;
	int biarc_flag = 0;
	double L1 = -10;
	double L2 = -10;

	double k1, k2;

	double L3, L4, E1, E2;

	double check = (psi - (theta2 - theta0) / 2);
	while (check > oct::math::pi || check < -1 * oct::math::pi)
	{
		if (check > oct::math::pi)
		{
			check = check - oct::math::pi;
		}
		else if (check < (-1 * oct::math::pi))
		{
			check = check + oct::math::pi;
		}
	}

	if (check > (-1 * smallAngle) && check < (smallAngle))
	{
		if ((psi - theta0) > (-1 * smallAngle) && (psi - theta0) < (smallAngle))
		{
			k1 = 0;
			k2 = k1;
			L1 = L / 2;
			L2 = L1;
		}
		else
		{
			k1 = -2 * sin((theta2 - theta0) / 2) / L;
			k2 = k1;
			L1 = std::abs(theta0 * L *
			              sin((oct::math::pi / 2 - theta0) / sin(2 * theta0)));
			L2 = L1;
		}
		estimateK = k1;
		estimatel = L1 + L2;

	} // if( check > (-1*smallAngle) && check < (smallAngle) )
	else
	{
		biarc_flag = 1;
		k1 = -4 * sin((3 * theta0 + theta2) / 4 - psi) *
		     cos((theta2 - theta0) / 4) / L;
		k2 = 4 * sin((theta0 + 3 * theta2) / 4 - psi) *
		     cos((theta2 - theta0) / 4) / L;
		double theta_join =
		    ComputeJoinTheta(x0, y0, theta0, x2, y2, theta2, k1, k2);

		if (k1 == 0)
		{
			L1 = L * sin((theta2 + theta0) / 2 - psi) /
			     sin((theta0 - theta2) / 2);
		}
		else
		{
			L1 = ComputeArcLength(theta0, theta_join, k1);
		}

		if (k2 == 0)
		{
			L2 = L * sin((theta2 + theta0) / 2 - psi) /
			     sin((theta0 - theta2) / 2);
		}
		else
		{
			L2 = ComputeArcLength(theta_join, theta2, k2);
		}

		double k3 = (4 / L) * cos((3 * theta0 + theta2) / 4 - psi) *
		            sin((theta2 - theta0) / 4);
		double k4 = (4 / L) * cos((theta0 + 3 * theta2) / 4 - psi) *
		            sin((theta2 - theta0) / 4);

		if (k3 == 0 && k4 == 0)
		{
			estimatel = L;
			estimateK = 0;
		}
		else
		{
			if (k3 != k4)
			{
				theta_join =
				    ComputeJoinTheta(x0, y0, theta0, x2, y2, theta2, k3, k4);

				if (k3 == 0)
				{
					L3 = L * sin((theta2 + theta0) / 2 - psi) /
					     sin((theta2 - theta0) / 2);
				}
				else
				{
					L3 = ComputeArcLength(theta0, theta_join, k3);
				}

				if (k4 == 0)
				{
					L4 = L * sin((theta2 + theta0) / 2 - psi) /
					     sin((theta0 - theta2) / 2);
				}
				else
				{
					L4 = ComputeArcLength(theta_join, theta2, k4);
				}
			}
		}

		E1 = (k2 - k1) * (k2 - k1);
		E2 = (k3 - k4) * (k3 - k4);

		if ((L1 + L2) < (L3 + L4))
		{
			estimateK = k1;
			estimatel = L1 + L2;
		}
		else
		{
			estimateK = k3;
			estimatel = L3 + L4;
		}
	}
	return 0;
}

double EulerSpiral::ComputeError(double x0,
                                           double y0,
                                           double theta0,
                                           double x2,
                                           double y2,
                                           double theta2,
                                           double k,
                                           double L)
{
	double ex, ey;
	double gamma = 2 * (theta2 - theta0 - k * L) / (L * L);
	EulerSpiralEndPoint(x0, y0, theta0, k, gamma, L, ex, ey);
	return sqrt((ex - x2) * (ex - x2) + (ey - y2) * (ey - y2));
}

int EulerSpiral::EulerSpiralEndPoint(double a,
                                               double b,
                                               double theta,
                                               double k,
                                               double gamma,
                                               double s,
                                               double & x,
                                               double & y)
{
	double epsilon = 0.00001;
	double constTerm = 0;
	double fc, fs, sc, ss;
	double cosTerm, sinTerm;
	double A, B;
	if ((gamma > 0 && gamma < epsilon) || (gamma < 0 && gamma > -1 * epsilon))
	{
		gamma = 0;
	}
	else
	{
		if (gamma > 0)
		{
			_private.fresnel((k + gamma * s) / (sqrt(oct::math::pi * gamma)), fc, fs);
			_private.fresnel((k) / (sqrt(oct::math::pi * gamma)), sc, ss);
			A = fc - sc;
			B = fs - ss;
			cosTerm = cos(theta - k * k * 0.5 / gamma);
			sinTerm = sin(theta - k * k * 0.5 / gamma);
			constTerm = sqrt(oct::math::pi / gamma);
			x = a + (constTerm) * (cosTerm * A - sinTerm * B);
			y = b + (constTerm) * (sinTerm * A + cosTerm * B);
		}
		if (gamma < 0)
		{
			_private.fresnel((k * -1 + gamma * s * -1) /
			            (sqrt(oct::math::pi * gamma * -1)),
			        fc,
			        fs);
			_private.fresnel(
			    (k * -1) / (sqrt(oct::math::pi * gamma * -1)), sc, ss);
			A = fc - sc;
			B = (fs - ss) * -1;
			cosTerm = cos(theta - k * k * 0.5 / gamma);
			sinTerm = sin(theta - k * k * 0.5 / gamma);
			constTerm = sqrt(oct::math::pi / gamma * -1);
			x = a + (constTerm) * (cosTerm * A - sinTerm * B);
			y = b + (constTerm) * (sinTerm * A + cosTerm * B);
		}
	}
	if (gamma == 0)
	{
		if (k == 0)
		{
			x = a + s * cos(theta);
			y = b + s * sin(theta);
		}
		else
		{
			constTerm = 0.1 / k;
			x = a + (constTerm) * (sin(k * s + theta) - sin(theta));
			y = b + (-1 * constTerm) * (cos(k * s + theta) - cos(theta));
		}
	}
	return 0;
}

int EulerSpiral::SolveIteratively(double x0,
                                            double y0,
                                            double theta0,
                                            double x1,
                                            double y1,
                                            double theta2,
                                            double estimateK,
                                            double estimateL,
                                            int iternum,
                                            double & Kfin,
                                            double & Lfin)
{
	double error =
	    ComputeError(x0, y0, theta0, x1, y1, theta2, estimateK, estimateL);
	double prevError = 1000;
	double errorStep = 0.1;
	double epsilon = 0.001;
	double epsilonError = 0.1;
	Kfin = estimateK;
	Lfin = estimateL;
	int i = 0;

	double error0, error1, error2, error3;
	error0 = error1 = error2 = error3 = 999;
	do
	{
		if (i == (iternum - 1) || error < epsilon)
		{
			return 1;
		}

		error0 = ComputeError(
		    x0, y0, theta0, x1, y1, theta2, Kfin + errorStep, Lfin);
		error1 = ComputeError(
		    x0, y0, theta0, x1, y1, theta2, Kfin - errorStep, Lfin);
		error2 = ComputeError(
		    x0, y0, theta0, x1, y1, theta2, Kfin, Lfin + errorStep);
		int f = 0;
		if (Lfin > errorStep)
		{
			f = 1;
			error3 = ComputeError(
			    x0, y0, theta0, x1, y1, theta2, Kfin, Lfin - errorStep);
		}
		if (f == 1)
		{
			error = _private.min4(error0, error1, error2, error3);
		}
		else
		{
			error = _private.min3(error0, error1, error2);
		}

		if (error >= prevError)
		{
			errorStep = errorStep / 2;
		}
		else if (error == error0)
		{
			Kfin = Kfin + errorStep;
		}
		else if (error == error1)
		{
			Kfin = Kfin - errorStep;
		}
		else if (error == error2)
		{
			Lfin = Lfin + errorStep;
		}
		else if (Lfin > errorStep)
		{
			Lfin = Lfin - errorStep;
		}
		prevError = error;
		i++;
		oct::serv::debug("DEBUG (Line ",
		                 __LINE__,
		                 "): Error : ",
		                 error0,
		                 " ",
		                 error1,
		                 " ",
		                 error2,
		                 " ",
		                 error3);
	} while (i <= iternum);

	return 0;
}

int EulerSpiral::SolveEuler(double x0,
                            double y0,
                            double theta0,
                            double x2,
                            double y2,
                            double theta2,
                            int iternum,
                            double & Kfin,
                            double & Lfin,
                            double & Gfin)
{
	double ek, el;
	int biarcRet =
	    _private.ComputeBiarcSolution(x0, y0, x2, y2, theta0, theta2, ek, el);

	int solveRet = SolveIteratively(
	    x0, y0, theta0, x2, y2, theta2, ek, el, iternum, Kfin, Lfin);
	Gfin = 2 * (theta2 - theta0 - Kfin * Lfin) / (Lfin * Lfin);

	return biarcRet || solveRet;
}

/* @Brief : Function to set the no. of divisions for numerical integraion of
 * fresnel integration This is required as for large difference b/w (x2,y2)
 * and (x0,y0) the default value 1000 is not enough
 * @params  [in] no_of_integral_steps    : No. of divisions for numerical
 * integraion of fresnel integration
 */

inline void EulerSpiral::setFresenlIntegrationSteps(int no_of_integral_steps)
{
	_private.integral_steps = no_of_integral_steps;
}
