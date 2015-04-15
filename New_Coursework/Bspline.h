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

	T differential(const float& t) const
	{
		T output;
		output = output * 0.0f;
		
		Pascal<ORDER> pascal;
		for (int i = 0; i < ORDER; ++i) {

			// differentiate using product rule
			float v1 = pow(t,i);
			float v2 = pow(1 - t, spline_order - i - 2);
			float v3 = i - spline_order + 1;
			float v4 = i * pow(t, i - 1);
			float v5 = pow(1 - t, spline_order - i - 1);

			// Two product protions, cut out zeroes as it breaks things.
			float second = ((i == 0.0f) ? 0.0f : v4*v5);
			float first = ((i == spline_order - 1) ? 0.0f : v1*v2*v3);

			// find value of this entry and add to the result
			T temp = control_points[i] * pascal(i)* (first + second);

			output += temp;
		}
		return output;
	}

	T dif(const float& t) const
	{
		T output;
		output = output * 0.0f;
		Pascal<ORDER> pascal;
		static const PascalTriangle<ORDER> pascl;
		for (size_t i = 0; i < ORDER; ++i)
		{
			int pow_single = i;
			int pow_brack = (ORDER - 1) - i;
			float u_mul = 0.0f;

			for (size_t j = 0; j <= pow_brack; ++j)
			{
				u_mul += (j + pow_single) * pascl(pow_brack)(j)*(j % 2 ? -1 : 1)*pow(t, j + pow_single - 1);
			}
			output += control_points[i] * u_mul;
		}
		return output;
	}

	inline const size_t order() const {
		return ORDER;
	}
private:
	T control_points[ORDER];
	static const int spline_order = ORDER;
};