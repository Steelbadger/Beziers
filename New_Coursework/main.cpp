#include<iostream>
#include<vector>
#include "Bspline.h"
#include "Vector3.h"
#include "Vector2.h"
#include "Bsurface.h"
#include "SplineInterface.h"

void main()
{
	Vector3 control_points[3] = { Vector3(0.0f,0.0f,0.0f), Vector3(0.0f, 0.5f, 0.0f), Vector3(0.0f, 1.0f, 0.25f) };
	Bspline<Vector3> my_spline(control_points);

	std::cout << "Vector3" << std::endl;
	std::cout << "Start: " << my_spline.at(0.0f) << std::endl;
	std::cout << "Half: " << my_spline.at(0.5f) << std::endl;
	std::cout << "End: " << my_spline.at(1.0f) << std::endl;
	std::cout << "Tangent at 0: " << my_spline.tangent(0.0f) << std::endl;
	std::cout << "Expected: " << Vector3(0.0f, 0.5f, 0.0f) - Vector3(0.0f, 0.0f, 0.0f) << std::endl << std::endl;


	Vector2 v2_points[4] = { Vector2(1.0f, 0.0f), Vector2(3.0f, 3.0f), Vector2(5.0f, 5.0f) , Vector2(7.0f, 2.0f)};
	Bspline<Vector2, 4> my_v2spline(v2_points);
	std::cout << "Vector2" << std::endl;
	std::cout << "Start: " << my_v2spline.at(0.0f) << std::endl;
	std::cout << "Quarter: " << my_v2spline.at(0.25f) << std::endl;
	std::cout << "End: " << my_v2spline.at(1.0f) << std::endl << std::endl;

	Vector3 tut1_s_1[3] = { Vector3(0.0f, 0.0f, 0.0f), Vector3(3,2,0), Vector3(4,0,0) };
	Vector3 tut1_s_2[3] = { Vector3(-1,4,2), Vector3(2,3.5f, 6), Vector3(4,3,2) };
	Vector3 tut1_s_3[3] = { Vector3(0,4,0), Vector3(2,3,0), Vector3(4,6,0) };

	Bspline<Vector3, 3> splines[3] = { tut1_s_1, tut1_s_2, tut1_s_3 };
	Bsurface<3> tut1_surf(splines);


	std::cout << "Surface" << std::endl;
	std::cout << "(0,0): " << tut1_surf.at(0, 0) << std::endl;
	std::cout << "(0,0.5): " << tut1_surf.at(0, 0.5) << std::endl;
	std::cout << "(0,1): " << tut1_surf.at(0, 1) << std::endl;
	std::cout << "(0.5,1): " << tut1_surf.at(0.5, 1) << std::endl;
	std::cout << "(1,1): " << tut1_surf.at(1, 1) << std::endl;
	std::cout << "(1,0.5): " << tut1_surf.at(1, 0.5) << std::endl;
	std::cout << "(1,0): " << tut1_surf.at(1, 0) << std::endl;
	std::cout << "(0.5,0): " << tut1_surf.at(0.5, 0) << std::endl;
	std::cout << "(0.5,0.5): " << tut1_surf.at(0.5, 0.5) << std::endl << std::endl;

	SplineInterface<Vector3>* temp = (SplineInterface<Vector3>*)&my_spline;
	std::cout << "Testing Temp" << std::endl;
	std::cout << "Start: " << temp->at(0.0f) << std::endl;
	std::cout << "Half: " << temp->at(0.5f) << std::endl;
	std::cout << "End: " << temp->at(1.0f) << std::endl << std::endl;

}

std::vector<Vector3> CreateTriangles(Bsurface<3> surf, size_t divisions)
{
	std::vector<Vector3> points;
	float precision = 1.0f / divisions;
	for (float u = 0; u <= 1.0f; u += precision)
	{
		for (float v = 0.0f; v <= 1.0f; v += precision)
		{
			points.push_back(surf.at(u, v));
		}
	}
	return points;
}