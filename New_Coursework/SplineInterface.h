#pragma once

template<class T>
class SplineInterface
{
public:
	virtual T at(const float& t) const = 0;
	virtual T tangent(const float& t) const = 0;
	virtual T dif(const float& t) const = 0;
	virtual const size_t order() const = 0;
};