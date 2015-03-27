#pragma once
#include "Bspline.h"
#include "Vector3.h"

template<size_t order = 3>
class Bsurface
{
public:
	Bsurface(Bspline<Vector3, order> controls[order])
	{
		for (size_t i = 0; i < order; ++i){
			control_splines[i] = controls[i];
		}
	}

	Bsurface(){};

	~Bsurface(){}

	Vector3 at(float u, float v)
	{
		Vector3 output(0.0f, 0.0f, 0.0f);
		Pascal<order> pascal;
		for (size_t i = 0; i < order; ++i)
		{
			int power1 = (order - 1) - i;
			int power2 = (order - 1) - (order - 1 - i);

			float one_mulp = pow(1 - u, power1);
			float sec_mulp = pow(u, power2);
			float pasc = pascal(i);

			output = output + control_splines[i].at(v) * pasc * one_mulp * sec_mulp;
		}
		return output;
	}

	Vector3 normal(float u, float v)
	{
		Vector3 tangent1 = { 0, 0, 0 };
		Vector3 tangent2 = { 0, 0, 0 };
		Vector3 norm;
		Pascal<order> pascals;
		for (size_t i = 0; i < order; ++i)
		{
			int power1 = (order - 1) - i;
			int power2 = (order - 1) - (order - 1 - i);
			float one_mulp = pow(1 - u, power1);
			float sec_mulp = pow(u, power2);
			float pasc = pascals(i);

			tangent1 += control_splines[i].tangent(v) * (pasc * one_mulp * sec_mulp);
			tangent2 += control_splines[i].at(u) * (pasc * pow(v, power1) * -power2 * pow(1 - v, power2 - 1) + pow(1 - v, power2) * pasc * power1 * pow(v, power1 - 1));
		}
		return Cross(tangent1, tangent2).Normalize();
	}

private:
	Bspline<Vector3, order> control_splines[order];
};