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
		static const PascalTriangle<order> pascl;
		for (size_t i = 0; i < order; ++i)
		{
			int power1 = (order - 1) - i;
			int power2 = i;
			float one_mulp = pow(1 - u, power1);
			float sec_mulp = pow(u, power2);
			float pasc = pascals(i);

			int pow_single = i;
			int pow_brack = (order - 1) - i;
			float u_mul= 0.0f;

			for (size_t j = 0; j <= pow_brack; ++j)
			{
				u_mul += (j + pow_single) * pascl(pow_brack)(j)*(j % 2 ? -1 : 1)*pow(u, j + pow_single - 1);
			}

			tangent1 += control_splines[i].at(v)*u_mul;
			tangent2 += control_splines[i].dif(v) * (pasc * one_mulp * sec_mulp);
		}
		return Cross(tangent1, tangent2).Normalize();
	}

private:
	Bspline<Vector3, order> control_splines[order];
};