#pragma once


template<size_t row>
class Pascal
{
public:
	Pascal()
	{
		if (!initialised) {
			initialised = true;
			for (size_t i = 1; i < row; ++i) {
				row_vals[i] = row_vals[i - 1] * ((row) - i) / float(i);
			}
		}
	}
	~Pascal(){}

	inline const unsigned int operator()(size_t element) const
	{
		return row_vals[element];
	}

private:
	static bool initialised;
	static unsigned int row_vals[row];
};

template<size_t row>
unsigned int Pascal<row>::row_vals[row] = { 1 };

template<size_t row>
bool Pascal<row>::initialised = false;