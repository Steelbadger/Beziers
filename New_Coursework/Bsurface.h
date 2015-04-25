#pragma once
#include "Bspline.h"
#include "Vector3.h"


class BSurfaceInterface
{
public:
	virtual Vector3 at(const float& u, const float&v) const = 0;
	virtual Vector3 normal(const float& u, const float& v) const = 0;
	virtual const unsigned int ord() const = 0;
};

template<size_t order = 3>
class Bsurface : public BSurfaceInterface
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

	Vector3 at(const float& u,const float& v) const
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

	Vector3 normal(const float& u,const float& v) const
	{
		Vector3 tangent1 = { 0, 0, 0 };
		Vector3 tangent2 = { 0, 0, 0 };
		Vector3 norm;
		Pascal<order> pascals;
		for (int i = 0; i < order; ++i)
		{
			int power1 = (order - 1) - i;
			int power2 = (order - 1) - (order - 1 - i);
			float pasc = pascals(i);

			float fir = (power2 == 0) ? 0.0f : pow(u, power1) * -power2 * pow(1 - u, power2 - 1);
			float sec = (power1 == 0) ? 0.0f : pow(1 - u, power2) * power1 * pow(u, power1 - 1);



			tangent1 += control_splines[i].differential(v) * pasc * (fir + sec);


			// differentiate using product rule
			float v1 = pow(u, i);
			float v2 = pow(1 - u, int(order) - i - 2);
			float v3 = i - int(order) + 1;
			float v4 = i * pow(u, i - 1);
			float v5 = pow(1 - u, int(order) - i - 1);

			// Two product protions, cut out zeroes as it breaks things.
			float second = ((i == 0.0f) ? 0.0f : v4*v5);
			float first = ((i == int(order) - 1) ? 0.0f : v1*v2*v3);

			// find value of this entry and add to the result
			tangent2 += control_splines[i].at(v) * pascals(i)* (first + second);
		}
		return Cross(tangent1, tangent2).Normalize();
	}

	const unsigned int ord() const
	{
		return order;
	}

private:
	Bspline<Vector3, order> control_splines[order];
};