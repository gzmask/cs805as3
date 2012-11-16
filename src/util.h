#ifndef UTIL_H
#define UTIL_H

//define preprocessing vars
#define IMG_COLS 512
#define IMG_ROWS 512
#define COLS 128 //X
#define ROWS 128 //Y
#define SLCS 128 //Z
#define VOL_LEN ( COLS * ROWS * SLCS )
#define IMG_LEN ( IMG_COLS * IMG_ROWS )

#include <array>
#include <vector>
#include <functional>
#include <iostream>
#include <string>
#include <iomanip>
#include <stdio.h>

//types
typedef std::array<unsigned char, IMG_LEN> ImagePanel;
typedef std::array<unsigned char, VOL_LEN> Volume;
typedef std::array<double, 3> Point;
typedef std::array<double, 2> Point2D;
typedef std::array<double, 3> Vector;
typedef std::array<int, 4> Four_counter;
typedef std::array<Vector, 3> UVN;
typedef std::vector<Point> Intersection;
/*
typedef struct {
	Point t0;
	Point t1;
} Intersection;
*/
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
void foreach_pixel_exec(ImagePanel*, std::function<unsigned char(Ray, Intersection, Volume*, Volume*)>, Volume*, Volume*);
void init_img_panel(ImagePanel*);
Ray ray_construction(int, int);
void compute_shading_volume(Volume*, Volume*);
bool ray_box_intersection(Ray, Intersection*);
unsigned char volume_ray_tracing(Ray, Intersection, Volume*, Volume*);
void compute_shading_volume(Volume*);
unsigned char tri_linear_interpolate(Volume*, Point);
unsigned char linear_interpolate(unsigned char, unsigned char, double);

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
std::array<int, 3> to_3d(int);
void print_img_panel(ImagePanel);
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
void save_to_file(ImagePanel*, std::string);
void print_ct_volume(Volume*);
void read_from_file(std::string, Volume*);
Point add(Point, Vector, double);

//global vars
extern Matrix Mwc;
extern Matrix Rwc;
extern Matrix Twc;
extern Matrix Mcw;
extern Matrix Rcw;
extern Matrix Tcw;
extern double xmin;
extern double ymin;
extern double xmax;
extern double ymax;
extern Point VRP;
extern Vector VPN;
extern Vector VUP;
extern double focal;
extern Vector Light;
extern double Ip;
#endif
