#pragma once

template<class T>
class SplineInterface
{
public:
	virtual T at(const float& t) = 0;
	virtual const size_t ord() = 0;
};