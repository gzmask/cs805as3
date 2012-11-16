#include "util.h"
#include <cmath>

//
bool ray_box_intersection(Ray ray, Intersection* intersection) {
  double t,x,y,z;
  int n = 0;

  //x=0
  t = -ray.ref[0]/ray.direction[0];
  y = ray.ref[1] + ray.direction[1]*t;
  z = ray.ref[2] + ray.direction[2]*t;
  if (y >0 && y<ROWS && z>0 && z<SLCS) {
    n++;
    intersection->push_back( {0, y, z} );
  }

  //x=127
  t = (127-ray.ref[0])/ray.direction[0];
  y = ray.ref[1] + ray.direction[1]*t;
  z = ray.ref[2] + ray.direction[2]*t;
  if (y >0 && y<ROWS && z>0 && z<SLCS) {
    n++;
    intersection->push_back( {127, y, z} );
  }
  
  //y=0
  t = -ray.ref[1]/ray.direction[1];
  x = ray.ref[0] + ray.direction[0]*t;
  z = ray.ref[2] + ray.direction[2]*t;
  if (x >0 && x<COLS && z>0 && z<SLCS) {
    n++;
    intersection->push_back( {x, 0, z} );
  }

  //y=127
  t = (127-ray.ref[1])/ray.direction[1];
  x = ray.ref[0] + ray.direction[0]*t;
  z = ray.ref[2] + ray.direction[2]*t;
  if (x >0 && x<COLS && z>0 && z<SLCS) {
    n++;
    intersection->push_back( {x, 127, z} );
  }

  //z=0
  t = -ray.ref[2]/ray.direction[2];
  x = ray.ref[0] + ray.direction[0]*t;
  y = ray.ref[1] + ray.direction[1]*t;
  if (x >0 && x<COLS && y>0 && y<ROWS) {
    n++;
    intersection->push_back( {x, y, 0} );
  }

  //z=127
  t = (127-ray.ref[2])/ray.direction[2];
  x = ray.ref[0] + ray.direction[0]*t;
  y = ray.ref[1] + ray.direction[1]*t;
  if (x >0 && x<COLS && y>0 && y<ROWS) {
    n++;
    intersection->push_back( {x, y, 127} );
  }

  if (n==2) {
    if (closer((*intersection)[1], (*intersection)[0], VRP)) {
      Point tp = (*intersection)[1];
      (*intersection)[1] = (*intersection)[0];
      (*intersection)[0] = tp;
    }
    return true;
  } else {
    return false;
  }
}


//the ray tracing function
unsigned char volume_ray_tracing(Ray ray, Intersection intersection, Volume* ct, Volume* color) {
  double Dt = 20.0; //the interval for sampling along the ray
  double C = 0.0; //for accumulating the shading value
  double T = 1.0; //for accumulating the transparency

  //marching
  for ( Point t = intersection[0];
      closer(t, intersection[1], intersection[0]);
      t = add(t, ray.direction, Dt) ) {
    //std::cout<<" x: "<<t[0]<<" y: "<<t[1]<<" z: "<<t[2]<<std::endl;
    auto a_i = (double)(tri_linear_interpolate(ct, t)/255.0);
    auto c_i = tri_linear_interpolate(color, t);
    C += (c_i * a_i);
    T *= (1.0 - c_i);
    if (T < 0)
      break;
  }

  //std::cout<<(int)C<<std::endl;
  return C;
}

//tri linear interpolation
unsigned char tri_linear_interpolate(Volume* vol, Point p) {
  int x = ceil(p[0]);
  int x_ = floor(p[0]);
  int y = ceil(p[1]);
  int y_ = floor(p[1]);
  int z = ceil(p[2]);
  int z_ = floor(p[2]);
  double px = p[0] - floor(p[0]);
  double py = p[1] - floor(p[1]);
  double pz = p[2] - floor(p[2]);

  return linear_interpolate(
    linear_interpolate(
      linear_interpolate((*vol)[to_1d(x,y,z)], (*vol)[to_1d(x_,y,z)], px),
      linear_interpolate((*vol)[to_1d(x,y_,z)], (*vol)[to_1d(x_,y_,z)], px),
      py),
    linear_interpolate(
      linear_interpolate((*vol)[to_1d(x,y,z_)], (*vol)[to_1d(x_,y,z_)], px),
      linear_interpolate((*vol)[to_1d(x,y_,z_)], (*vol)[to_1d(x_,y_,z_)], px),
      py),
    pz);
}

//linear interploation
unsigned char linear_interpolate(unsigned char a, unsigned char b, double p) {
  auto diff = b - a;
  return (unsigned char)(a+diff*p);
}

//shading the ct volume
void compute_shading_volume(Volume* ct, Volume* color) {
  for (int i=0; i<VOL_LEN; i++) {
    std::array<int, 3> p = to_3d(i); 
    //ignore the outest layer
    if (p[0] == 0 || p[0] == (COLS-1) ||
        p[1] == 0 || p[1] == (COLS-1) ||
        p[2] == 0 || p[2] == (COLS-1)) {
      (*color)[i] = 0;
      continue;
    }

    Vector N_ = {
      (double)((*ct)[i] - (*ct)[to_1d(p[0]-1, p[1], p[2])]), 
      (double)((*ct)[i] - (*ct)[to_1d(p[0], p[1]-1, p[2])]), 
      (double)((*ct)[i] - (*ct)[to_1d(p[0], p[1], p[2]-1)])
    };

    //threshold the blur surfaces
    if (get_length(N_) < 10) {
      (*color)[i] = 0;
      continue;
    }

    //calculate shading value
    Vector N = normalize(N_);
    unsigned char I = Ip * dot_product(N, Light);
    (*color)[i] = I;
  }
  return;
}

//pixel iterator for img panel.
void foreach_pixel_exec(ImagePanel* img, std::function<unsigned char(Ray, Intersection, Volume*, Volume*)> ray_func, Volume* ct, Volume* color) {
  for (int i = 0; i < IMG_LEN; i++) { //for each pixel
    //using to_2d function to get x,y camera coordinates
    std::array<int, 2> cam_xy = to_2d(i);

    //construct Ray
    Ray ray = ray_construction(cam_xy[0], cam_xy[1]);

    //get intersection
    Intersection* intersection = new Intersection;
    if (ray_box_intersection(ray, intersection)) {
      //get shading value using the passed-in functor
      (*img)[i] = ray_func(ray, *intersection, ct, color);
    }
    delete intersection;
  }
  return;
}

//ray constructor
Ray ray_construction(int x, int y) {
  //calculate x unit
  double x_delta = (xmax-xmin) / IMG_COLS;
  double y_delta = (ymax-ymin) / IMG_ROWS;

  //calculate the point on img panel with world coordinate
  double x_ = xmax - x_delta * x;
  double y_ = ymax - y_delta * y;

  //get vector v0. it is trival that VRP is p0
  Point p0 = VRP;
  Point p1_ = {x_, y_, focal};
  Point p1 = mul(Mcw, p1_);
  Vector v0_ = {p1[0] - p0[0],
    p1[1] - p0[1], 
    p1[2] - p0[2]};
  Vector v0 = normalize(v0_);

  return { p0[0], p0[1], p0[2],
           v0[0], v0[1], v0[2]};
}

//initialize img panel to all 0s
void init_img_panel(ImagePanel* img) {
  for (auto& pixel: *img) { //foreach pixel in empty_img
    pixel = 0;
  }
  return;
}

//==========helper functions==========

//stepping a point on some direction
Point add(Point t, Vector N, double Dt) {
  Vector N_ = {N[0]*Dt, N[1]*Dt, N[2]*Dt};
  return {t[0]+N_[0], t[1]+N_[1], t[2]+N_[2]};
}

//read file to array
void read_from_file(std::string filename, Volume* ct) {
  FILE *infid;

  //open file
  if ((infid = fopen(filename.c_str(), "rb")) == NULL) {
    std::cout<<"Open CT DATA File Error."<<std::endl;
    return; 
  }
  std::cout<<"Sucessfully open CT DATA File."<<std::endl;

  //read file
  for (int i=0; i<SLCS; i++) { /* Read one slice at a time. */
    int n = fread(&(*ct)[to_1d(0,0,i)], sizeof(char), COLS*ROWS, infid);
    if (n < COLS*ROWS*sizeof(char)) {
      std::cout<<"Read CT data slice "<<i<<" error."<<std::endl;
      return; 
    }
  }
  std::cout<<"Sucessfully read CT DATA File."<<std::endl;
  fclose(infid);
  return; 
}

//save image panel to binary file
void save_to_file(ImagePanel* img, std::string filename) {
  FILE *pFile;
  pFile = fopen ( filename.c_str(), "wb" );
  std::cout<<"img size written to output.raw: "<<img->size()<<std::endl;
  fwrite (img->data(), sizeof(int), img->size() , pFile );
  fclose (pFile);
  return;
}

//get the D of Ax+By+Cz+D=0 from POLY4
double get_D_poly4(POLY4 obj) {
  return -(obj.N[0]*obj.v1[0]+obj.N[1]*obj.v1[1]+obj.N[2]*obj.v1[2]);
}

//check if point is inside a 4 side polygon
bool in_poly4(Point p, POLY4 obj) {
  //flatten the polygon and the point
  POLY4_2D obj2d = flatten(obj, p);

  //tranlate flatten polygon to origin
  obj2d.v1[0] = obj2d.v1[0] - obj2d.p[0];
  obj2d.v1[1] = obj2d.v1[1] - obj2d.p[1];
  obj2d.v2[0] = obj2d.v2[0] - obj2d.p[0];
  obj2d.v2[1] = obj2d.v2[1] - obj2d.p[1];
  obj2d.v3[0] = obj2d.v3[0] - obj2d.p[0];
  obj2d.v3[1] = obj2d.v3[1] - obj2d.p[1];
  obj2d.v4[0] = obj2d.v4[0] - obj2d.p[0];
  obj2d.v4[1] = obj2d.v4[1] - obj2d.p[1];

  //count intersections with v=0 (u>0) (u<0) and u=0 (v>0) (v<0)
  Four_counter counter = {0,0,0,0};
  counter = count_intersection(obj2d.v1, obj2d.v2, counter);
  counter = count_intersection(obj2d.v2, obj2d.v3, counter);
  counter = count_intersection(obj2d.v3, obj2d.v4, counter);
  counter = count_intersection(obj2d.v4, obj2d.v1, counter);

  //std::cout<<"count: "<<counter[0]<<", "<<counter[1]<<", "<<counter[2]<<", "<<counter[3]<<std::endl;

  if (counter[0] %2 != 0 && counter[1] %2 != 0 && counter[2] %2 != 0 && counter[3] %2 != 0)
    return true;
  else
    return false;
}

//count intersections of edge v1-v2 with u and v axises
Four_counter count_intersection(Point2D v1, Point2D v2, Four_counter counter) {
  int u_plus_count = counter[0];
  int u_minus_count = counter[1];
  int v_plus_count = counter[2];
  int v_minus_count = counter[3];

  //if v2[1]-v1[1] is 0
  if (v2[1]-v1[1] == 0) {
    /*           
             |
             |v1--v2
        -----|---->
             |
             |
    */
    if (v2[0]*v1[0] > 0) {
      return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
    }
    /*           
             |
         v1--|--v2
        -----|---->
             |
             |
    */
    if (v2[0]*v1[0] < 0) {
      if (v1[1] >= 0) {
        v_plus_count++;
        return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
      } else {
        v_minus_count++;
        return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
      }
    }
  }

  //if v2[0]-v1[0] is 0
  if (v2[0]-v1[0] == 0) {
    /*           
         v1 |
         |  |
         v2 |
        ----|--->
            |
            |
    */
    if (v2[1]*v1[1] > 0) {
      return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
    }
    /*           
         v1 |
         |  |
        -|--|--->
         |  |
         v2 |
    */
    if (v2[1]*v1[1] < 0) {
      if (v1[0] >= 0) {
        u_plus_count++;
        return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
      } else {
        u_minus_count++;
        return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
      }
    }
  }

  //first calcualte the slope m = (y2 − y1)/(x2 − x1)
  double m = (v2[1] - v1[1]) / (v2[0] - v1[0]);

  //when v = 0, u = -y1/m + x1
  /*           
       v1     |
        \     |
      ---*----|--->
          \   |
           v2 |
  */
  double u = - v1[1]/m + v1[0];
  Point2D intersection1 = {u, 0};
  if (inside_bounding(v1, v2, intersection1)) {
    if (u>=0) {
      return {u_plus_count+1, u_minus_count, v_plus_count, v_minus_count};
    } else {
      return {u_plus_count, u_minus_count+1, v_plus_count, v_minus_count};
    }
  }

  //when u = 0, v = -x1*m + y1
  /*           
         |
     ----|--->
       v1|
        \|
         *
         |\
         | v2  
  */
  double v = - v1[0]*m + v1[1];
  Point2D intersection2 = {0, v};
  if (inside_bounding(v1, v2, intersection2)) {
    if (v>=0) {
      return {u_plus_count, u_minus_count, v_plus_count+1, v_minus_count};
    } else {
      return {u_plus_count, u_minus_count, v_plus_count, v_minus_count+1};
    }
  }

  return {u_plus_count, u_minus_count, v_plus_count, v_minus_count};
}

//check if third point is inside the bounding box from point 1 and 2
bool inside_bounding(Point2D v1, Point2D v2, Point2D p) {
  double u_max, u_min, v_max, v_min;

  if (v1[0] >= v2[0]) {
    u_max = v1[0];
    u_min = v1[0];
  } else {
    u_max = v2[0];
    u_min = v2[0];
  }

  if (v1[1] >= v2[1]) {
    v_max = v1[1];
    v_min = v1[1];
  } else {
    v_max = v2[1];
    v_min = v2[1];
  }

  if ( (p[0] > u_min) && (p[0] < u_max) && 
       (p[1] > v_min) && (p[1] < v_max) )
    return true;
  else
    return false;
}

//make 2D polygon from 3D polygon
POLY4_2D flatten(POLY4 obj, Point p) {
  //find out the dominated dimension
  int drop_i = find_max(obj.N[0], obj.N[1], obj.N[2]);

  //drop the dominated dimension
  if (drop_i == 0) {
    POLY4_2D result = {
      obj.v1[1], obj.v1[2],
      obj.v2[1], obj.v2[2],
      obj.v3[1], obj.v3[2],
      obj.v4[1], obj.v4[2],
      p[1], p[2]
      };
    return result;
  } else if (drop_i == 1) {
    POLY4_2D result = {
      obj.v1[0], obj.v1[2],
      obj.v2[0], obj.v2[2],
      obj.v3[0], obj.v3[2],
      obj.v4[0], obj.v4[2],
      p[0], p[2]
      };
    return result;
  } else if (drop_i == 2) {
    POLY4_2D result = {
      obj.v1[0], obj.v1[1],
      obj.v2[0], obj.v2[1],
      obj.v3[0], obj.v3[1],
      obj.v4[0], obj.v4[1],
      p[0], p[1]
      };
    return result;
  }

  POLY4_2D null_ = { -1,-1, -1,-1, -1,-1, -1,-1 };
  return null_;
}

//find out the max index for three doubles, 0: first, 1: second, 2: third
int find_max(double x, double y, double z) {
  if (fabs(x) >= fabs(y) && fabs(x) >= fabs(z)) {
    return 0;
  } else if (fabs(y) >= fabs(x) && fabs(y) >= fabs(z)) {
    return 1;
  } else if (fabs(z) >= fabs(x) && fabs(z) >= fabs(y)) {
    return 2;
  } else {
    return -1;
  }
}

//get transformation matrix
Matrix get_T(Point vrp) {
  Row r1 = {1, 0, 0, -vrp[0]};
  Row r2 = {0, 1, 0, -vrp[1]};
  Row r3 = {0, 0, 1, -vrp[2]};
  Row r4 = {0, 0, 0, 1};
  return {r1, r2, r3, r4};
}

//get inverse transformation matrix
Matrix get_Ti(Point vrp) {
  Row r1 = {1, 0, 0, vrp[0]};
  Row r2 = {0, 1, 0, vrp[1]};
  Row r3 = {0, 0, 1, vrp[2]};
  Row r4 = {0, 0, 0, 1};
  return {r1, r2, r3, r4};
}

//get rotation matrix
Matrix get_R(Point vrp, Vector vpn, Vector vup) {
  //first get the translation matrix from world to view
  //auto mt = get_T(vrp);

  //we can see vpn_ and vup_ as vectors. such that we can apply them to get_uvn function from q2
  auto uvn = get_uvn(vpn, vup);
  //finally contruct our roation matrix using method 2 on class notes
  Row r1 = { uvn[0][0],uvn[0][1],uvn[0][2],0 };
  Row r2 = { uvn[1][0],uvn[1][1],uvn[1][2],0 };
  Row r3 = { uvn[2][0],uvn[2][1],uvn[2][2],0 };
  Row r4 = { 0, 0, 0, 1 };
  return { r1, r2, r3, r4 };
}

//get inverse rotation matrix
Matrix get_Ri(Point vrp, Vector vpn, Vector vup) {
  Matrix m = get_R(vrp, vpn, vup);
  Row r1 = { m[0][0], m[1][0], m[2][0], m[3][0] };
  Row r2 = { m[0][1], m[1][1], m[2][1], m[3][1] };
  Row r3 = { m[0][2], m[1][2], m[2][2], m[3][2] };
  Row r4 = { m[0][3], m[1][3], m[2][3], m[3][3] };
  return {r1,r2,r3,r4};
}

//world to camera
Matrix get_M(Point vrp, Vector vpn, Vector vup) {
  return mul(get_R(vrp, vpn, vup), get_T(vrp));
}

//camera to world
Matrix get_Mi(Point vrp, Vector vpn, Vector vup) {
  return mul(get_Ti(vrp), get_Ri(vrp, vpn, vup));
}


//matrix multiplication
Point mul(Matrix m, Point x) {
  return mul(x, m);
}

Point mul(Point x, Matrix m) {
  double w =  m[3][0] * x[0]
        + m[3][1] * x[1]
        + m[3][2] * x[2]
        + m[3][3];
  return {(x[0]*m[0][0]+x[1]*m[0][1]+x[2]*m[0][2]+m[0][3])/w,
          (x[0]*m[1][0]+x[1]*m[1][1]+x[2]*m[1][2]+m[1][3])/w,
          (x[0]*m[2][0]+x[1]*m[2][1]+x[2]*m[2][2]+m[2][3])/w};
}

Matrix mul(Matrix m, Matrix n) {
  Row r1 = {m[0][0]*n[0][0]+m[0][1]*n[1][0]+m[0][2]*n[2][0]+m[0][3]*n[3][0],
            m[0][0]*n[0][1]+m[0][1]*n[1][1]+m[0][2]*n[2][1]+m[0][3]*n[3][1],
            m[0][0]*n[0][2]+m[0][1]*n[1][2]+m[0][2]*n[2][2]+m[0][3]*n[3][2],
            m[0][0]*n[0][3]+m[0][1]*n[1][3]+m[0][2]*n[2][3]+m[0][3]*n[3][3]};
  Row r2 = {m[1][0]*n[0][0]+m[1][1]*n[1][0]+m[1][2]*n[2][0]+m[1][3]*n[3][0],
            m[1][0]*n[0][1]+m[1][1]*n[1][1]+m[1][2]*n[2][1]+m[1][3]*n[3][1],
            m[1][0]*n[0][2]+m[1][1]*n[1][2]+m[1][2]*n[2][2]+m[1][3]*n[3][2],
            m[1][0]*n[0][3]+m[1][1]*n[1][3]+m[1][2]*n[2][3]+m[1][3]*n[3][3]};
  Row r3 = {m[2][0]*n[0][0]+m[2][1]*n[1][0]+m[2][2]*n[2][0]+m[2][3]*n[3][0],
            m[2][0]*n[0][1]+m[2][1]*n[1][1]+m[2][2]*n[2][1]+m[2][3]*n[3][1],
            m[2][0]*n[0][2]+m[2][1]*n[1][2]+m[2][2]*n[2][2]+m[2][3]*n[3][2],
            m[2][0]*n[0][3]+m[2][1]*n[1][3]+m[2][2]*n[2][3]+m[2][3]*n[3][3]};
  Row r4 = {m[3][0]*n[0][0]+m[3][1]*n[1][0]+m[3][2]*n[2][0]+m[3][3]*n[3][0],
            m[3][0]*n[0][1]+m[3][1]*n[1][1]+m[3][2]*n[2][1]+m[3][3]*n[3][1],
            m[3][0]*n[0][2]+m[3][1]*n[1][2]+m[3][2]*n[2][2]+m[3][3]*n[3][2],
            m[3][0]*n[0][3]+m[3][1]*n[1][3]+m[3][2]*n[2][3]+m[3][3]*n[3][3]};
  return {r1,r2,r3,r4};
}

Row mul(Row x, Matrix m) {
  return {x[0]*m[0][0]+x[1]*m[0][1]+x[2]*m[0][2]+x[3]*m[0][3],
          x[0]*m[1][0]+x[1]*m[1][1]+x[2]*m[1][2]+x[3]*m[1][3],
          x[0]*m[2][0]+x[1]*m[2][1]+x[2]*m[2][2]+x[3]*m[2][3],
          x[0]*m[3][0]+x[1]*m[3][1]+x[2]*m[3][2]+x[3]*m[3][3]};
}

Row mul(Matrix m, Row x) {
  return mul(x, m);
}

//return if p1 is closer to p0 than p2
bool closer(Point p1, Point p2, Point p0) {
  Vector v1 = { (p1[0] - p0[0]), (p1[1] - p0[1]), (p1[2] - p0[2])};
  double d1 = get_length(v1);
  Vector v2 = { (p2[0] - p0[0]), (p2[1] - p0[1]), (p2[2] - p0[2])};
  double d2 = get_length(v2);
  return d1 < d2;
}

//Translate 2D array index of row column to 1D index.
//Notice that x, or column index, starts with 0. 
//If return value is -1 then there is an out-of-bounce error.
int to_1d(int x, int y) {
  if (x >= IMG_COLS || x < 0)
    return -1;
  if (y >= IMG_ROWS || y < 0)
    return -1;
  return (IMG_COLS*y + x);
}

//Translate 3D array index to 1D index.
//x is counting in column, y is counting in row, z is counting in slice
int to_1d(int x, int y, int z) {
  if (x >= COLS || x < 0)
    return -1;
  if (y >= ROWS || y < 0)
    return -1;
  if (z >= SLCS || z < 0)
    return -1;
  return (COLS*ROWS*z + COLS*y + x);
}

//Translate 1d array index to 2d
std::array<int, 2> to_2d(int x) {
  if (x>=(IMG_COLS*IMG_ROWS) || x < 0) {
    return {-1,-1};
  }
  int y_ = x / IMG_COLS; 
  int x_ = x % IMG_COLS;
  return {x_, y_};
}

std::array<int, 3> to_3d(int x) {
  if (x>=(COLS*ROWS*SLCS) || x < 0) {
    return {-1,-1, -1};
  }
  int z_ = x / (COLS*ROWS);
  int y_ = (x % (COLS*ROWS)) / COLS;
  int x_ = x % COLS;
  return {x_, y_, z_};
}

//prints the ct volume
void print_ct_volume(Volume* ct) {
  std::cout<<std::endl;
  for (auto& voxel : *ct) {
    std::cout<<(int)voxel<<", ";
  }
  std::cout<<std::endl<<"Volume array size: "<<ct->size()<<std::endl;
}

//prints the img panel
void print_img_panel(ImagePanel img) {
  std::cout<<std::endl;
  for (auto& pixel : img) {
    std::cout<<(int)pixel<<", ";
  }
  std::cout<<std::endl<<"Array size: "<<img.size()<<std::endl;
}

//get u,v,n from two non-collinear vectors
UVN get_uvn(Vector V1, Vector V2) {

  //get n, which is just normalized V1
  Vector n = normalize(V1); 

  //get u, which is normalized V2 x V1
  Vector u = normalize(cross_product(V2, V1));

  //get v, which is normalized n x u
  Vector v = normalize(cross_product(n, u));

  return {u,v,n};
}

//normalize a Vector
Vector normalize(Vector x) {
  return { x[0]/get_length(x), 
           x[1]/get_length(x), 
           x[2]/get_length(x) }; 
}

//dot product
double dot_product(Vector x, Vector y) {
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; 

}

//calculates cross product of two Vectors
Vector cross_product(Vector x, Vector y) {
  return { x[1]*y[2] - x[2]*y[1],
           x[2]*y[0] - x[0]*y[2],
           x[0]*y[1] - x[1]*y[0]};
}

//calculates length of a Vector
double get_length(Vector x) {
  return sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
}

