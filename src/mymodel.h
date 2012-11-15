#include <stdio.h>

/* Definition of image buffers */

#define SLCS 128
#define ROWS 128
#define COLS 128
unsigned char	CT[SLCS][ROWS][COLS]; /* a 3D array for CT data */
unsigned char	SHADING[SLCS][ROWS][COLS]; /* a 3D array for shading values */

#define IMG_ROWS 512
#define IMG_COLS 512
unsigned char	out_img[IMG_ROWS][IMG_COLS];

/* Camera parameters */
float VRP[3] = {128.0, 64.0, 250.0};
float VPN[3] = {-64.0, 0.0, -186.0};
float VUP[3] = {0.0, 1.0, 0.0};

/* Image Plane Sizes */
float focal = 0.05;	/* 50 mm lens */
float xmin = -0.0175;	/* 35 mm "film" */
float ymin = -0.0175;
float xmax = 0.0175;
float ymax = 0.0175;

/* Light direction (unit length vector) */
float Light[3] = {0.577, -0.577, -0.577};
/* Light Intensity */
float Ip = 255.0;
