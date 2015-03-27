#pragma once
#include<vector>


namespace meta
{
	// defines the condition in which the class diverges to base case definition.
	#define PASCAL_CONDITION(pArg0, pArg1)	\
		(pArg1) == 0 || (pArg0) == (pArg1)

	// the pascal identity value.
	#define PASCAL_IDENTITY(pArg0, pArg1)									\
		pascal_identity<PASCAL_CONDITION(pArg0, pArg1), pArg0, pArg1>::value

	// the pascal identity: (n; k) = (n - 1; k - 1) + (n - 1; k)
	template <bool B, size_t N, size_t K>
	struct pascal_identity
	{	
		// defines a class-wide constant value for a particular template instance.
		// the definition is recursive, so other template instances are generated until the
		// base case is reached. the constant values are then added at compile time.
		enum {value = PASCAL_IDENTITY(N - 1, K - 1) + PASCAL_IDENTITY(N - 1, K)};
	};

	//explicit class definition - base case: bool = true.
	template <size_t N, size_t K>
	struct pascal_identity<true, N, K>
	{
		// class-wide constant value is always 1 whenever the first template argument is true.
		enum {value = 1};
	};

	template<size_t row, size_t element>
	struct Pascal
	{
		enum {value = PASCAL_IDENTITY(row, element)};
	};

#undef PASCAL_CONDITION
#undef PASCAL_IDENTITY
};


class LamePascal
{
public:
	LamePascal(){}
	LamePascal(size_t row)
	{
		vals.resize(row, 1);
		for (size_t i = 1; i < row; ++i) {
			vals[i] = vals[i - 1] * ((row)-i) / float(i);
		}
	}

	~LamePascal(){};


	inline const unsigned int operator()(size_t element) const
	{
		return vals[element];
	}
private:
	std::vector<unsigned int> vals;
};


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


template<size_t rows>
class PascalTriangle
{
public:
	PascalTriangle(){
		for (size_t i = 0; i < rows; ++i) {
			all_rows[i] = LamePascal(i + 1);
		}
	}
	~PascalTriangle() {}

	inline const LamePascal& operator()(size_t element) const
	{
		return all_rows[element];
	}

private:
	LamePascal all_rows[rows];
};