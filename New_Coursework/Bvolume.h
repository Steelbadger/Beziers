#pragma once
#include "Bsurface.h"
#include "Vector3.h"


class BVolumeInterface
{
public:
	virtual Vector3 at(const float& u, const float&v, const float& s) const = 0;
	virtual const unsigned int ord() const = 0;
};

template<size_t order = 3>
class Bvolume : public BVolumeInterface
{
public:
	Bvolume(Bsurface<order> controls[order])
	{
		for (size_t i = 0; i < order; ++i){
			control_surfaces[i] = controls[i];
		}
	}

	Bvolume(){};

	~Bvolume(){}

	Vector3 at(const float& u, const float& v, const float& s) const
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

			output = output + control_surfaces[i].at(v,s) * pasc * one_mulp * sec_mulp;
		}
		return output;
	}

	const unsigned int ord() const
	{
		return order;
	}

private:
	Bsurface<order> control_surfaces[order];
};