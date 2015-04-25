#include<iostream>
#include<vector>
#include "Bspline.h"
#include "Vector3.h"
#include "Vector2.h"
#include "Bsurface.h"
#include "SplineInterface.h"
#include "Bvolume.h"

#include "glut\glut.h"
#include <stdlib.h>
#include <stdio.h>
#include <Windows.h>


#ifndef CALLBACK
#define CALLBACK
#endif

const int ORDER = 4;

Vector3 ctlpoints[ORDER][ORDER];
Vector3 volpoints[ORDER][ORDER][ORDER];

enum Axis {X_AXIS, Y_AXIS, Z_AXIS};

Axis influence = Z_AXIS;
int pSel = 0;
int sSel = 0;
int fSel = 0;
float rotation = 0.0f;
int showPoints = 0;
float d = 3.0f;
float r = d / 2.0f;
std::vector<Vector3> tris;
std::vector<Vector3> norm;

std::vector<Vector3> meshverts;
std::vector<Vector3> meshnorm;

/*
*  Initializes the control points of the surface to a small hill.
*  The control points range from -3 to +3 in x, y, and z
*/

bool LoadObj(const char*);
std::vector<Vector3> Spherize(std::vector<Vector3>);
std::vector<Vector3> SubDivide(std::vector<Vector3> m);
void UnitCube();

void init_surface(void)
{
	int u, v;
	for (u = 0; u < ORDER; u++) {
		for (v = 0; v < ORDER; v++) {
			float z;

			int difu = u - (ORDER / 2);
			int difv = v - (ORDER / 2);

			z = 3.0f - ((float(difu) / ORDER) * (float(difv) / ORDER) * (3.0f*4.0f));

			if ((u == 1 || u == 2) && (v == 1 || v == 2))
				z = 3.0;
			else
				z = -3.0;

			ctlpoints[u][v] = Vector3(2.0f*(((float(u) / ORDER)*4.0f) - 1.5), 2.0f*(((float(v) / ORDER)*4.0f) - 1.5), z);

		}
	}
}

void init_volume()
{
	int u, v, z;

	for (u = 0; u < ORDER; u++) {
		for (v = 0; v < ORDER; v++) {
			for (z = 0; z < ORDER; z++) {
				volpoints[u][v][z] = Vector3((((float(u) / (ORDER-1))*d) - r), (((float(v) / (ORDER-1))*d) - r), (((float(z) / (ORDER-1))*d) - r));
			}
		}
	}
}

void CreateTriangles(const BSurfaceInterface* surf, size_t divisions)
{
	std::vector<Vector3> points;
	std::vector<Vector3> normals;
	float precision = 1.0f / divisions;
	for (float u = 0; u <= 1.0f; u += precision)
	{
		for (float v = 0.0f; v <= 1.0f; v += precision)
		{
			points.push_back(surf->at(u, v));
			normals.push_back(surf->normal(u, v));
		}
	}

	std::vector<Vector3> triangles;
	std::vector<Vector3> norms;

	for (int i = 0; i < divisions ; i++) {
		for (int j = 0; j < divisions ; j++) {
			triangles.push_back(points[i + j*divisions]);
			triangles.push_back(points[1 + i + (1 + j)*divisions]);
			triangles.push_back(points[i + (1 + j)*divisions]);


			triangles.push_back(points[1 + i + (1 + j)*divisions]);
			triangles.push_back(points[i + j*divisions]);
			triangles.push_back(points[1 + i + j*divisions]);

			Vector3 f1 = Normal(points[i + j*divisions], points[1 + i + (1 + j)*divisions], points[i + (1 + j)*divisions])*-1;
			Vector3 f2 = Normal(points[1 + i + (1 + j)*divisions], points[i + j*divisions], points[1 + i + j*divisions])*-1;

			norms.push_back(f1);
			norms.push_back(f1);
			norms.push_back(f1);
			norms.push_back(f2);
			norms.push_back(f2);
			norms.push_back(f2);
		}
	}

	tris = triangles;
	norm = norms;
}

void CreateTriangles(const BVolumeInterface* surf, std::vector<Vector3> fpoints)
{
	std::vector<Vector3> points;
	for (auto i = fpoints.begin(); i != fpoints.end(); ++i)
	{
		points.push_back(surf->at(i->x+0.5, i->y+0.5, i->z+0.5));
	}

	std::vector<Vector3> norms;

	for (int i = 0; i < points.size(); i += 3){
		Vector3 f1 = Normal(points[i], points[i+1], points[i+2]);
		norms.push_back(f1);
		norms.push_back(f1);
		norms.push_back(f1);
	}

	tris = points;
	norm = norms;
}

void RecalculateSurface()
{
	Bspline<Vector3, ORDER> norm_splines[ORDER];

	for (int i = 0; i < ORDER; ++i){
		norm_splines[i] = ctlpoints[i];
	}

	Bsurface<ORDER> norm_surf(norm_splines);

	CreateTriangles(&norm_surf, 100);
}

void RecalculateVolume()
{
	Bspline<Vector3, ORDER> norm_splines[ORDER][ORDER];

	for (int i = 0; i < ORDER; ++i){
		for (int j = 0; j < ORDER; ++j){
			norm_splines[i][j] = volpoints[i][j];
		}
	}

	Bsurface<ORDER> norm_surf[ORDER];
	for (int i = 0; i < ORDER; ++i){

			norm_surf[i] = norm_splines[i];

	}
	
	Bvolume<ORDER> norm_vol(norm_surf);

	CreateTriangles(&norm_vol, meshverts);
}

void TestSplines()
{
	Bspline<Vector3, ORDER> norm_splines[ORDER];

	for (int i = 0; i < ORDER; ++i){
		norm_splines[i] = ctlpoints[i];
	}

	Bsurface<ORDER> norm_surf(norm_splines);

	std::cout << "result at 0.5, 0.5: " << norm_surf.at(0.5f, 0.5f)<<std::endl;
}

void TestVolume()
{
	Bspline<Vector3, ORDER> norm_splines[ORDER][ORDER];

	for (int i = 0; i < ORDER; ++i){
		for (int j = 0; j < ORDER; ++j){
			norm_splines[i][j] = volpoints[i][j];
		}
	}

	Bsurface<ORDER> norm_surf[ORDER];
	for (int i = 0; i < ORDER; ++i){

		norm_surf[i] = norm_splines[i];

	}

	Bvolume<ORDER> norm_vol(norm_surf);

	std::cout << "result at 0.0, 0.0, 0.0: " << norm_vol.at(0.0f, 0.0f, 0.0f) << std::endl;
	std::cout << "result at 1.0, 0.0, 0.0: " << norm_vol.at(1.0f, 0.0f, 0.0f) << std::endl;
	std::cout << "result at 0.0, 1.0, 0.0: " << norm_vol.at(0.0f, 1.0f, 0.0f) << std::endl;
	std::cout << "result at 0.0, 0.0, 1.0: " << norm_vol.at(0.0f, 0.0f, 1.0f) << std::endl;
}


/*  Initialize material property and depth buffer.
*/
void init(void)
{
	GLfloat mat_diffuse[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 100.0 };

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glFrontFace(GL_CCW);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_CULL_FACE);

	glCullFace(GL_BACK);

//	init_surface();
	
//	RecalculateSurface();
//	TestSplines();

	UnitCube();
	meshverts = SubDivide(meshverts);
	meshverts = SubDivide(meshverts);
	meshverts = SubDivide(meshverts);
	meshverts = SubDivide(meshverts);
	meshverts = Spherize(meshverts);
	meshnorm = meshverts;
	init_volume();
	RecalculateVolume();
	TestVolume();
}

void display(void)
{
	GLfloat knots[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
	glRotatef(330.0f, 1.0f, 0.0f, 0.0f);
	glRotatef(rotation, 0.0f, 0.0f, 1.0f);
	glScalef(0.25, 0.25, 0.25);

	glEnable(GL_LIGHTING);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < tris.size(); ++i) {
		glVertex3fv((GLfloat*)&tris[i]);
		glNormal3fv((GLfloat*)&norm[i]);
	}
	glEnd();



	if (showPoints) {
		glPointSize(5.0);
		glDisable(GL_LIGHTING);
		glColor3f(1.0, 1.0, 0.0);
		glScalef(2.0f, 2.0f, 2.0f);
		glBegin(GL_POINTS);
		for (int i = 0; i < ORDER; i++) {
			for (int j = 0; j < ORDER; j++) {
				for (int k = 0; k < ORDER; k++){
					glVertex3fv((GLfloat*)&volpoints[i][j][k]);
				}
			}
		}
		glEnd();

		glPointSize(8.0);
		glColor3f(1.0, 0.0, 0.0);
		glDisable(GL_DEPTH_TEST);
		glBegin(GL_POINTS);
		glVertex3fv((GLfloat*)&volpoints[fSel][sSel][pSel]);
		glEnd();
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
	}
	glPopMatrix();
	glFlush();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (GLdouble)w / (GLdouble)h, 3.0, 8.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0, 0.0, -5.0);
}

void update()
{

}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'c':
	case 'C':
		showPoints = !showPoints;
		glutPostRedisplay();
		break;
	case 'x':
	case 'X':
		influence = X_AXIS;
		std::cout << "Moving on X-Axis" << std::endl;
		break;
	case 'y':
	case 'Y':
		influence = Y_AXIS;
		std::cout << "Moving on Y-Axis" << std::endl;
		break;
	case 'z':
	case 'Z':
		influence = Z_AXIS;
		std::cout << "Moving on Z-Axis" << std::endl;
		break;
	case 'a':
	case 'A':
		// increment selected point
		pSel = (++pSel) % ORDER;
		glutPostRedisplay();
		break;
	case 'd':
	case 'D':
		// decrement selected point
		pSel = (--pSel) % ORDER;
		pSel = pSel < 0 ? ORDER + pSel : pSel;
		glutPostRedisplay();
		break;
	case 'e':
	case 'E':
		// increment selected point
		fSel = (++fSel) % ORDER;
		glutPostRedisplay();
		break;
	case 'q':
	case 'Q':
		// decrement selected point
		fSel = (--fSel) % ORDER;
		fSel = fSel < 0 ? ORDER + fSel : fSel;
		glutPostRedisplay();
		break;
	case 'w':
	case 'W':
		// increment spline
		sSel = (++sSel) % ORDER;
		glutPostRedisplay();
		break;
	case 's':
	case 'S':
		// decrement spline
		sSel = (--sSel) % ORDER;
		sSel = sSel < 0 ? ORDER + sSel : sSel;
		glutPostRedisplay();
		break;
	case 'r':
	case 'R':
		// randomise position
		{
			float x = (((float(fSel) / (ORDER - 1))*d) - r) + ((float(rand() % 600) - 300.0f) / 300.0f)*r / ORDER;
			float y = (((float(sSel) / (ORDER - 1))*d) - r) + ((float(rand() % 600) - 300.0f) / 300.0f)*r / ORDER;
			float z = (((float(pSel) / (ORDER - 1))*d) - r) + ((float(rand() % 600) - 300.0f) / 300.0f)*r / ORDER;

			volpoints[fSel][sSel][pSel] = Vector3(x, y, z);
		}
		RecalculateVolume();
		glutPostRedisplay();
		break;
	case 't':
	case 'T':
		for (int i = 0; i < ORDER; ++i){
			for (int j = 0; j < ORDER; ++j){
				for (int k = 0; k < ORDER; ++k) {
					float x = (((float(i) / (ORDER - 1))*d) - r) + ((float(rand() % 600) - 300.0f) / 300.0f)*r / ORDER;
					float y = (((float(j) / (ORDER - 1))*d) - r) + ((float(rand() % 600) - 300.0f) / 300.0f)*r / ORDER;
					float z = (((float(k) / (ORDER - 1))*d) - r) + ((float(rand() % 600) - 300.0f) / 300.0f)*r / ORDER;

					volpoints[i][j][k] = Vector3(x,y,z);
				}
			}
		}
//		RecalculateSurface();
		RecalculateVolume();
		glutPostRedisplay();
		break;
	case 'b':
	case 'B':
		init_volume();
//		RecalculateSurface();
		RecalculateVolume();
		glutPostRedisplay();
		break;
	case ',':
		rotation += 5.0f;
		glutPostRedisplay();
		break;
	case '.':
		rotation -= 5.0f;
		glutPostRedisplay();
		break;
	case ']':
		if (influence == X_AXIS)
			volpoints[fSel][sSel][pSel].x += 0.2f;
		if (influence == Y_AXIS)
			volpoints[fSel][sSel][pSel].y += 0.2f;
		if (influence == Z_AXIS)
			volpoints[fSel][sSel][pSel].z += 0.2f;

//		RecalculateSurface();
		RecalculateVolume();
		glutPostRedisplay();
		break;
	case '[':
		if (influence == X_AXIS)
			volpoints[fSel][sSel][pSel].x -= 0.2f;
		if (influence == Y_AXIS)
			volpoints[fSel][sSel][pSel].y -= 0.2f;
		if (influence == Z_AXIS)
			volpoints[fSel][sSel][pSel].z -= 0.2f;
//		RecalculateSurface();
		RecalculateVolume();
		glutPostRedisplay();
		break;

	case 27:
		exit(0);
		break;
	default:
		break;
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	init();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMainLoop();
	return 0;
}


//void main()
//{
//	Vector3 control_points[3] = { Vector3(0.0f,0.0f,0.0f), Vector3(0.0f, 0.5f, 0.0f), Vector3(0.0f, 1.0f, 0.25f) };
//	Bspline<Vector3> my_spline(control_points);
//
//	std::cout << "Vector3" << std::endl;
//	std::cout << "Start: " << my_spline.at(0.0f) << std::endl;
//	std::cout << "0.25: " << my_spline.at(0.25f) << std::endl;
//	std::cout << "End: " << my_spline.at(1.0f) << std::endl;
//	std::cout << "Differential at 0.1: " << my_spline.differential(0.1f) << std::endl;
//	std::cout << "Expected: " << Vector3(0.0f, 0.5f, 0.0f) - Vector3(0.0f, 0.0f, 0.0f) << std::endl << std::endl;
//
//
//	Vector2 v2_points[4] = { Vector2(1.0f, 0.0f), Vector2(3.0f, 3.0f), Vector2(5.0f, 5.0f) , Vector2(7.0f, 2.0f)};
//	Bspline<Vector2, 4> my_v2spline(v2_points);
//	std::cout << "Vector2" << std::endl;
//	std::cout << "Start: " << my_v2spline.at(0.0f) << std::endl;
//	std::cout << "Quarter: " << my_v2spline.at(0.25f) << std::endl;
//	std::cout << "End: " << my_v2spline.at(1.0f) << std::endl << std::endl;
//
//	Vector3 tut1_s_1[3] = { Vector3(0.0f, 0.0f, 0.0f), Vector3(3,2,0), Vector3(4,0,0) };
//	Vector3 tut1_s_2[3] = { Vector3(-1,4,2), Vector3(2,3.5f, 6), Vector3(4,3,2) };
//	Vector3 tut1_s_3[3] = { Vector3(0,4,0), Vector3(2,3,0), Vector3(4,6,0) };
//
//	Bspline<Vector3, 3> splines[3] = { tut1_s_1, tut1_s_2, tut1_s_3 };
//	Bsurface<3> tut1_surf(splines);
//
//	Pascal<1> first;
//	Pascal<2> second;
//	Pascal<3> third;
//	Pascal<4> fourth;
//	first(0);
//	std::cout << "\t\t\t" << first(0) << std::endl;
//	std::cout << "\t\t" << second(0) << "\t\t" << second(1) << std::endl;
//	std::cout << "\t" << third(0) << "\t\t" << third(1) << "\t\t" << third(2) << std::endl;
//	std::cout << "" << fourth(0) << "\t\t" << fourth(1) << "\t\t" << fourth(2) << "\t\t" << fourth(3) << std::endl;
//
//	std::cout << "Surface" << std::endl;
//	std::cout << "(0,0): " << tut1_surf.at(0, 0) << std::endl;
//	std::cout << "(0,0.5): " << tut1_surf.at(0, 0.5) << std::endl;
//	std::cout << "(0,1): " << tut1_surf.at(0, 1) << std::endl;
//	std::cout << "(0.5,1): " << tut1_surf.at(0.5, 1) << std::endl;
//	std::cout << "(1,1): " << tut1_surf.at(1, 1) << std::endl;
//	std::cout << "(1,0.5): " << tut1_surf.at(1, 0.5) << std::endl;
//	std::cout << "(1,0): " << tut1_surf.at(1, 0) << std::endl;
//	std::cout << "(0.5,0): " << tut1_surf.at(0.5, 0) << std::endl;
//	std::cout << "(0.5,0.5): " << tut1_surf.at(0.5, 0.5) << std::endl << std::endl;
//
//	std::cout << "Normals:" << std::endl;
//	std::cout << "(0, 0): " << tut1_surf.normal(0.2, 0.3) << std::endl << std::endl;
//
//	SplineInterface<Vector3>* temp = (SplineInterface<Vector3>*)&my_spline;
//	std::cout << "Testing Temp" << std::endl;
//	std::cout << "Start: " << temp->at(0.0f) << std::endl;
//	std::cout << "Half: " << temp->at(0.5f) << std::endl;
//	std::cout << "End: " << temp->at(1.0f) << std::endl << std::endl;
//
//
//	Vector3 norm_1[3] = { Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(2, 0, 0) };
//	Vector3 norm_2[3] = { Vector3(0, 1, 0), Vector3(1, 1, 3), Vector3(2, 1, 0) };
//	Vector3 norm_3[3] = { Vector3(0, 2, 0), Vector3(1, 2, 0), Vector3(2, 2, 0) };
//
//	Bspline<Vector3, 3> norm_splines[3] = { norm_1, norm_2, norm_3 };
//	Bsurface<3> norm_surf(norm_splines);
//
//	std::cout << "Test Normals: " << std::endl;
//	std::cout << "(0.43, 0.88)" << norm_surf.normal(0.43f, 0.88f) << std::endl;
//
//}



std::vector<Vector3> CreateNormals(Bsurface<3> surf, size_t divisions)
{
	std::vector<Vector3> points;
	float precision = 1.0f / divisions;
	for (float u = 0; u <= 1.0f; u += precision)
	{
		for (float v = 0.0f; v <= 1.0f; v += precision)
		{
			points.push_back(surf.normal(u, v));
		}
	}
	return points;
}

bool LoadObj(const char* path)
{
	std::vector<unsigned int> vertIndices, uvIndices, normalIndices;
	std::vector<Vector3> tempVerts;
	std::vector<Vector3> tempNormals;
	std::vector<Vector2> tempUVs;

	FILE * file;
	fopen_s(&file, path, "r");
	if (file == NULL){
		std::cout << "Cannot Open File: " << path << std::endl;
		return false;
	}

	float max = 0;
	float min = 0;

	while (true){

		char lineHeader[128];
		// read the first word of the line
		int res = fscanf_s(file, "%s", lineHeader);
		if (res == EOF)
			break; // EOF = End Of File. Quit the loop.

		// else : parse lineHeader

		if (strcmp(lineHeader, "v") == 0){
			Vector3 vertex;
			fscanf_s(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
			tempVerts.push_back(vertex);
			max = std::fmaxf(vertex.x, max);
			max = std::fmaxf(vertex.y, max);
			max = std::fmaxf(vertex.z, max);

			min = std::fminf(vertex.x, min);
			min = std::fminf(vertex.y, min);
			min = std::fminf(vertex.z, min);
		}
		else if (strcmp(lineHeader, "vt") == 0){
			Vector2 uv;
			fscanf_s(file, "%f %f\n", &uv.u, &uv.v);
			tempUVs.push_back(uv);
		}
		else if (strcmp(lineHeader, "vn") == 0){
			Vector3 normal;
			fscanf_s(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
			tempNormals.push_back(normal);
		}
		else if (strcmp(lineHeader, "f") == 0){
			std::string vertex1, vertex2, vertex3;
			unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
			int matches = fscanf_s(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2]);
			if (matches != 9){
				std::cout << "Cannot Read File: " << path << std::endl;
				fclose(file);
				return false;
			}
			vertIndices.push_back(vertexIndex[0]);
			vertIndices.push_back(vertexIndex[1]);
			vertIndices.push_back(vertexIndex[2]);
			uvIndices.push_back(uvIndex[0]);
			uvIndices.push_back(uvIndex[1]);
			uvIndices.push_back(uvIndex[2]);
			normalIndices.push_back(normalIndex[0]);
			normalIndices.push_back(normalIndex[1]);
			normalIndices.push_back(normalIndex[2]);
		}
		else {
			// Probably a comment, eat up the rest of the line
			char stuff[1000];
			fgets(stuff, 1000, file);
		}

	}

	// For each vertex of each triangle
	for (unsigned int i = 0; i<vertIndices.size(); i++){

		// Get the indices of its attributes
		unsigned int vertexIndex = vertIndices[i];
		unsigned int uvIndex = uvIndices[i];
		unsigned int normalIndex = normalIndices[i];

		// Get the attributes thanks to the index
		Vector3 vertex = (tempVerts[vertexIndex - 1] - min)/(max-min);
		Vector2 uv = tempUVs[uvIndex - 1];
		Vector3 normal = tempNormals[normalIndex - 1];

		// Put the attributes in buffers
		meshverts.push_back(vertex);
	}
	fclose(file);
//	numVerts = vertIndices.size();

	return true;
}

std::vector<Vector3> SubDivide(std::vector<Vector3> m)
{
	//  Now create a bunch of new faces (non indexed)
	std::vector<Vector3> newPoints;

	for (int i = 0; i < m.size(); i += 3) {
		Vector3 p1 = m[i];
		Vector3 p2 = m[i + 1];
		Vector3 p3 = m[i + 2];
		Vector3 p12 = (p1 + p2) / 2;
		Vector3 p13 = (p1 + p3) / 2;
		Vector3 p23 = (p2 + p3) / 2;

		//  Push four new triangles to the newPoints array

		newPoints.push_back(p1);
		newPoints.push_back(p12);
		newPoints.push_back(p13);

		newPoints.push_back(p2);
		newPoints.push_back(p23);
		newPoints.push_back(p12);

		newPoints.push_back(p12);
		newPoints.push_back(p23);
		newPoints.push_back(p13);

		newPoints.push_back(p13);
		newPoints.push_back(p23);
		newPoints.push_back(p3);
	}

	return newPoints;
}

std::vector<Vector3> Spherize(std::vector<Vector3> m)
{
	std::vector<Vector3> out;
	for (auto i = m.begin(); i != m.end(); ++i)
	{
		out.push_back(i->Normalize());
	}
	return out;
}

void UnitCube()
{
	std::vector<Vector3> verts;


	std::vector<Vector3> points;

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				points.push_back(Vector3(-0.5 + i, -0.5 + j, -0.5 + k));
			}
		}
	}


	//  Face 1 - Triangle 1
	verts.push_back(points[1]);
	verts.push_back(points[0]);
	verts.push_back(points[5]);




	//  Face 1 - Triangle 2
	verts.push_back(points[5]);
	verts.push_back(points[0]);
	verts.push_back(points[4]);



	//  Face 2 - Triangle 1
	verts.push_back(points[5]);
	verts.push_back(points[4]);
	verts.push_back(points[7]);



	//  Face 2 - Triangle 2
	verts.push_back(points[7]);
	verts.push_back(points[4]);
	verts.push_back(points[6]);


	//  Face 3 - Triangle 1
	verts.push_back(points[3]);
	verts.push_back(points[2]);
	verts.push_back(points[1]);



	//  Face 3 - Triangle 2
	verts.push_back(points[1]);
	verts.push_back(points[2]);
	verts.push_back(points[0]);


	//  Face 4 - Triangle 1
	verts.push_back(points[0]);
	verts.push_back(points[2]);
	verts.push_back(points[4]);


	//  Face 4 - Triangle 2
	verts.push_back(points[4]);
	verts.push_back(points[2]);
	verts.push_back(points[6]);



	//  Face 5 - Triangle 1
	verts.push_back(points[2]);
	verts.push_back(points[3]);
	verts.push_back(points[6]);



	//  Face 5 - Triangle 2
	verts.push_back(points[6]);
	verts.push_back(points[3]);
	verts.push_back(points[7]);



	//  Face 6 - Triangle 1
	verts.push_back(points[3]);
	verts.push_back(points[1]);
	verts.push_back(points[7]);


	//  Face 6 - Triangle 2
	verts.push_back(points[7]);
	verts.push_back(points[1]);
	verts.push_back(points[5]);

	meshverts = verts;
}