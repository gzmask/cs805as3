/* IMPORTANT: After you download this file, you should rename it
	to "mymodel.h"
*/
typedef std::array<float, 4> Row;
typedef std::array<Row, 4> Matrix;

/* Definition of the structure for Sphere */
typedef struct {
	float x, y, z;	/* center of the circle */
	float radius;	/* radius of the circle */
	float kd;	/* diffuse reflection coefficient */
} SPHERE;

/* Definition of Polygon with 4 edges */
typedef struct {
	float v[4][3];	/* list of vertices */
	float N[3];	/* normal of the polygon */
	float kd;	/* diffuse reflection coefficient */
} POLY4;

/* create a spherical object */
SPHERE obj1 = {1.0, 1.0, 1.0,	/* center of the circle */
		 1.0,		/* radius of the circle */
		 0.75};		/* diffuse reflection coefficient */

/* create a polygon object */
POLY4 obj2 = {	0.0, 0.0, 0.0,	/* v0 */
		0.0, 0.0, 2.0,	/* v1 */
		2.0, 0.0, 2.0,	/* v2 */
		2.0, 0.0, 0.0,	/* v3 */
		0.0, 1.0, 0.0,	/* normal of the polygon */
		0.8};		/* diffuse reflection coefficient */

/* definition of the image buffer */
#define ROWS 512
#define COLS 512

unsigned char img[ROWS][COLS];

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */
float xmin = 0.0175;
float ymin = -0.0175;
float xmax = -0.0175;
float ymax = 0.0175;


/* definition of the camera parameters */
float VRP[3] = {1.0, 2.0, 3.5};
float VPN[3] = {0.0, -1.0, -2.5};
float VUP[3] = {0.0, 1.0, 0.0};

float focal = 0.05;	/* focal length simulating 50 mm lens */


/* definition of light source */
float LPR[3] = {-10.0, 10.0, 2.0};	/* light position */
float Ip = 200.0;	/* intensity of the point light source */

/* === transformation matrices (to be constructed) === */

/* Transformation from the world to the camera coordinates */
Matrix Mwc =
	{1.0, 0.0, 0.0, 0.0,
	 0.0, 1.0, 0.0, 0.0,
	 0.0, 0.0, 1.0, 0.0,
	 0.0, 0.0, 0.0, 1.0};

/* Transformation from the camera to the world coordinates */
Matrix Mcw =
	{1.0, 0.0, 0.0, 0.0,
	 0.0, 1.0, 0.0, 0.0,
	 0.0, 0.0, 1.0, 0.0,
	 0.0, 0.0, 0.0, 1.0};
