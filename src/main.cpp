#include <iostream>
#include "util.h"

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

//unsigned char img[ROWS][COLS];

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */
double xmin = 0.0175;
double ymin = -0.0175;
double xmax = -0.0175;
double ymax = 0.0175;

/* definition of the camera parameters */
//Point VRP = {1.0, 2.0, 3.5};
//Vector VPN = {0.0, -1.0, -2.5};
//Vector VUP = {0.0, 1.0, 0.0};
Point VRP = {1.0, 1.0, 4.5};
Vector VPN = {0.0, 0.0, -1.0};
Vector VUP = {0.0, 1.0, 0.0};

double focal = 0.05;	/* focal length simulating 50 mm lens */

/* definition of light source */
Point LRP = {-10.0, 10.0, 2.0};	/* light position */
double Ip = 200.0;	/* intensity of the point light source */

/* Transformation from the world to the camera coordinates */
Matrix Mwc = get_M(VRP, VPN, VUP);
Matrix Rwc = get_R(VRP, VPN, VUP);
Matrix Twc = get_T(VRP);
/* Transformation from the camera to the world coordinates */
Matrix Mcw = get_Mi(VRP, VPN, VUP);
Matrix Rcw = get_Ri(VRP, VPN, VUP);
Matrix Tcw = get_Ti(VRP);
/* Transformation from the world to light coordinates */
Matrix Mwl = get_T(LRP);
/* Transformation from the light to the world coordinates */
Matrix Mlw = get_Ti(LRP);

void unit_test() {
  //tests
  Point vrp = {6.0, 10.0, -5.0};
  Vector vpn = {-6.0, -9.0, 5.0}; 
  Vector vup = {0.0, 1.0, 0.0};
  auto uvn = get_uvn(vpn, vup);
  std::cout<<"get_uvn function:"<<std::endl;
  for (auto vecotr : uvn) {//for each Vecotr in uvn
    for (auto num : vecotr) {//for each number in Vecotr
      std::cout<<num<<',';
    }   
    std::cout<<std::endl;
  }
  
  Matrix mwc = get_M(vrp, vpn, vup);
  Matrix mcw = get_Mi(vrp, vpn, vup);
  Matrix twc = get_T(vrp);
  Matrix tcw = get_Ti(vrp);
  Matrix rwc = get_R(vrp, vpn, vup);
  Matrix rcw = get_Ri(vrp, vpn, vup);
  pmatrix("mwc:", mwc);
  pmatrix("twc:", twc);
  pmatrix("rwc:", rwc);
  pmatrix("mcw:", mcw);
  pmatrix("tcw:", tcw);
  pmatrix("rcw:", rcw);
  pmatrix("Mcw:", Mcw);

  std::cout<<"to_1d function, expected to be 512:"<<std::endl;
  std::cout<<to_1d(0, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be 0, 1:"<<std::endl;
  std::cout<<to_2d(512)[0]<<std::endl;
  std::cout<<to_2d(512)[1]<<std::endl;

  std::cout<<"to_1d function, expected to be 513:"<<std::endl;
  std::cout<<to_1d(1, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be 1,1:"<<std::endl;
  std::cout<<to_2d(513)[0]<<std::endl;
  std::cout<<to_2d(513)[1]<<std::endl;

  std::cout<<"to_1d function, expected to be 1023:"<<std::endl;
  std::cout<<to_1d(511, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be 511,1:"<<std::endl;
  std::cout<<to_2d(1023)[0]<<std::endl;
  std::cout<<to_2d(1023)[1]<<std::endl;

  std::cout<<"to_1d function, expected to be -1:"<<std::endl;
  std::cout<<to_1d(512, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be -1,-1:"<<std::endl;
  std::cout<<to_2d(512*512)[0]<<std::endl;
  std::cout<<to_2d(512*512)[1]<<std::endl;

  std::cout<<"closer function, expected to be 1 and 0:"<<std::endl;
  std::cout<<closer({1,1,1},{2,2,2},{0,0,0})<<std::endl;
  std::cout<<closer({3,3,3},{2,2,2},{0,0,0})<<std::endl;

  std::cout<<"mul function: point * matrix"<<std::endl;
  Point a = {3,3,3};
  std::cout<<mul(a, Mcw)[0]<<std::endl;
  std::cout<<mul(Mcw, a)[1]<<std::endl;
  std::cout<<mul(a, Mcw)[2]<<std::endl;

  pmatrix("mul funciton: matrix*matrix:", mul(Mwc, Mcw));

  std::cout<<"mul function: row * matrix"<<std::endl;
  Row b = {3,3,3,3};
  std::cout<<mul(b, Mcw)[0]<<std::endl;
  std::cout<<mul(Mcw, b)[1]<<std::endl;
  std::cout<<mul(b, Mcw)[2]<<std::endl;
  std::cout<<mul(b, Mcw)[3]<<std::endl;

  std::cout<<"to_1d, expected to be 165130:"<<to_1d(10, 10, 10)<<std::endl;
  std::cout<<"to_1d, expected to be -1:"<<to_1d(200, 10, 10)<<std::endl;
  std::cout<<"to_1d, expected to be -1:"<<to_1d(10, 200, 10)<<std::endl;
  std::cout<<"to_1d, expected to be -1:"<<to_1d(10, 10, 200)<<std::endl;
  std::cout<<"to_3d, expected to be 10,10,10:"<<to_3d(165130)[0]<<to_3d(165130)[1]<<to_3d(165130)[2]<<std::endl;
  std::cout<<"to_3d, expected to be 127,127,127:"<<to_3d(2097151)[0]<<to_3d(2097151)[1]<<to_3d(2097151)[2]<<std::endl;
  std::cout<<"to_3d, expected to be -1,-1,-1:"<<to_3d(2097153)[0]<<to_3d(2097153)[1]<<to_3d(2097153)[2]<<std::endl;
  std::cout<<"to_3d, expected to be -1,-1,-1:"<<to_3d(2097153)[0]<<to_3d(2097153)[1]<<to_3d(2097153)[2]<<std::endl;
  std::cout<<"volumn length: "<<VOL_LEN<<std::endl;

  std::cout<<"read_from_file function: "<<std::endl;
  CTVolume *tmpct = new CTVolume;
  read_from_file("smallHead.den", tmpct);
  //print_ct_volume(tmpct);
  delete tmpct;

  return;
}

int main () {

  unit_test();

  //main program for ray tracing
  /*
  ImagePanel img;
  img = init_img_panel(img);
  img = foreach_pixel_exec(img, ray_tracing);
  //print_img_panel(img);
  save_to_file(img, "output.raw"); //save result to binary file
  */

  //main program for volume rendering
  CTVolume* ct = new CTVolume;
  read_from_file("smallHead.den", ct);
  print_ct_volume(ct);
  compute_shading_volume(ct);
  print_ct_volume(ct);
  ImagePanel* img = new ImagePanel;
  init_img_panel(img);
  foreach_pixel_exec(img, volume_ray_tracing);
  save_to_file(img, "output.raw"); //save result to binary file
  delete img;
  delete ct;

  return 0;
}
