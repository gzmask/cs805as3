#ifndef UTIL_H
#define UTIL_H

//define preprocessing vars
#define IMG_X 512
#define IMG_Y 512
#define IMG_LEN ( IMG_X * IMG_Y )
/* definition of the image buffer */
#define ROWS IMG_Y
#define COLS IMG_X

#include <array>
#include <functional>
#include <iostream>
#include <string>
#include <iomanip>
#include <stdio.h>

//types
typedef std::array<unsigned char, IMG_LEN> ImagePanel;
typedef std::array<double, 3> Point;
typedef std::array<double, 2> Point2D;
typedef std::array<double, 3> Vector;
typedef std::array<int, 4> Four_counter;
typedef std::array<Vector, 3> UVN;
typedef struct {
	Point intersection;	/* intersection point */
	Vector normal;	/* intersection polygon normal vector */
	double kd;	/* diffuse reflection coefficient of the surface */
} Intersection;
typedef struct {
	Point ref;	/* reference point, where the ray is from */
	Vector direction;	/* ray direction */
} Ray;
typedef std::array<double, 4> Row;
typedef std::array<Row, 4> Matrix;
typedef struct {
	double x, y, z;	/* center of the circle */
	double radius;	/* radius of the circle */
	double kd;	/* diffuse reflection coefficient */
} SPHERE;
typedef struct {
        Point v1;       /* list of vertices */
        Point v2;
        Point v3;
        Point v4;
	Vector N;	/* normal of the polygon */
	double kd;	/* diffuse reflection coefficient */
} POLY4;
typedef struct {
        Point2D v1;
        Point2D v2;
        Point2D v3;
        Point2D v4;
        Point2D p;
} POLY4_2D;

//functions
POLY4_2D flatten(POLY4, Point);
ImagePanel foreach_pixel_exec(ImagePanel, std::function<unsigned char(Ray, Point)>);
ImagePanel init_img_panel(ImagePanel);
unsigned char ray_tracing(Ray, Point);
Intersection ray_objects_intersection(Ray);
unsigned char shading(Intersection, Point);
Intersection ray_sphere_intersection(Ray, SPHERE);
Intersection ray_polygon_intersection(Ray, POLY4);
Ray ray_construction(int, int);

//helper functions
bool inside_bounding(Point2D, Point2D, Point2D);
Four_counter count_intersection(Point2D, Point2D, Four_counter);
int find_max(double, double, double);
double get_D_poly4(POLY4);
bool in_poly4(Point, POLY4);
Point mul(Point, Matrix);
Point mul(Matrix, Point);
Matrix mul(Matrix, Matrix);
Row mul(Row, Matrix);
Row mul(Matrix, Row);
int to_1d(int, int);
int to_1d(int, int, int);
std::array<int, 2> to_2d(int);
void print_img_panel(ImagePanel);
void pmatrix(std::string, Matrix);
bool closer(Point, Point, Point);
UVN get_uvn(Vector V1, Vector V2);
Matrix get_T(Point);
Matrix get_Ti(Point);
Matrix get_R(Point, Vector, Vector);
Matrix get_Ri(Point, Vector, Vector);
Matrix get_M(Point, Vector, Vector);
Matrix get_Mi(Point, Vector, Vector);
double get_length(Vector);
Vector cross_product(Vector, Vector);
double dot_product(Vector, Vector);
Vector normalize(Vector);
void save_to_file(ImagePanel);
void read_from_file(std::string);

//global vars
extern Matrix Mwc;
extern Matrix Rwc;
extern Matrix Twc;
extern Matrix Mcw;
extern Matrix Rcw;
extern Matrix Tcw;
extern Matrix Mwl;
extern Matrix Mlw;
extern double xmin;
extern double ymin;
extern double xmax;
extern double ymax;
extern Point VRP;
extern Vector VPN;
extern Vector VUP;
extern double focal;
extern Point LRP;
extern double Ip;
extern SPHERE obj1;
extern POLY4 obj2;
#endif
