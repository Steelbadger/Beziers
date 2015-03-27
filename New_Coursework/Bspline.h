#pragma once
#include <math.h>
#include "Pascal.h"
#include "SplineInterface.h"

template<class T, size_t order = 3>
class Bspline: public SplineInterface<T>
{
public:
	Bspline(const T controls[order])
	{
		for (size_t i = 0; i < order; ++i){
			control_points[i] = controls[i];
		}
	};
	Bspline(){}
	 
	~Bspline(){};

	T at(const float& t)
	{
		T output;
		output = output * 0.0f;
		Pascal<order> pascal;
		for (size_t i = 0; i < order; ++i)
		{
			int power1 = (order - 1) - i;
			int power2 = (order - 1) - (order - 1 - i);

			float one_mulp = pow(1 - t, power1);
			float sec_mulp = pow(t, power2);
			float pasc = pascal(i);

			output = output + control_points[i] * pasc * one_mulp * sec_mulp;
		}
		return output;
	}

	T tangent(const float& t)
	{
		T output;
		output = output * 0.0f;
		Pascal<order> pascal;
		for (size_t i = 0; i < order; ++i)
		{
			int power1 = (order - 1) - i;
			int power2 = (order - 1) - (order - 1 - i);
			float pasc = pascal(i);

			output += control_points[i] * (pasc * pow(t, power1) * -power2 * pow(1 - t, power2 - 1) + pow(1 - t, power2) * pasc * power1 * pow(t, power1 - 1));
		}
		return output;
	}

	const size_t ord(){
		return spline_order;
	}
private:
	T control_points[order];
	static const size_t spline_order = order;
};