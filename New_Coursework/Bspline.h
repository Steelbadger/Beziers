#pragma once
#include <math.h>
#include "Pascal.h"
#include "SplineInterface.h"

template<class T, size_t ORDER = 3>
class Bspline: public SplineInterface<T>
{
public:
	Bspline(const T controls[ORDER])
	{
		for (size_t i = 0; i < ORDER; ++i){
			control_points[i] = controls[i];
		}
	};
	Bspline(){}
	 
	~Bspline(){};

	T at(const float& t) const
	{
		T output;
		output = output * 0.0f;
		Pascal<ORDER> pascal;
		for (size_t i = 0; i < ORDER; ++i)
		{
			int power1 = (ORDER - 1) - i;
			int power2 = (ORDER - 1) - (ORDER - 1 - i);

			float one_mulp = pow(1 - t, power1);
			float sec_mulp = pow(t, power2);
			float pasc = pascal(i);

			output = output + control_points[i] * pasc * one_mulp * sec_mulp;
		}
		return output;
	}

	T tangent(const float& t) const
	{
		T output;
		output = output * 0.0f;
		Pascal<ORDER> pascal;
		for (size_t i = 0; i < ORDER; ++i)
		{
			int power1 = (ORDER - 1) - i;
			int power2 = (ORDER - 1) - (ORDER - 1 - i);
			float pasc = pascal(i);

			output += control_points[i] * (pasc * pow(t, power1) * -power2 * pow(1 - t, power2 - 1) + pow(1 - t, power2) * pasc * power1 * pow(t, power1 - 1));
		}
		return output;
	}

	inline const size_t order() const {
		return ORDER;
	}
private:
	T control_points[ORDER];
	static const size_t spline_order = ORDER;
};