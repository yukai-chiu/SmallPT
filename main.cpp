#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<iostream>
#include<random>

std:: default_random_engine generator;
std::uniform_real_distribution<double> distr(0.0,1.0);

double erand48(int X){
    return distr(generator);
}

enum Refl_t {DIFF, SPEC, REFR};
class Vector3{
    public:
        double x, y, z;
        Vector3(double x_=0, double y_=0, double z_=0){
            x = x_;
            y = y_;
            z = z_;
        }
        Vector3 operator+(const Vector3 &b) const {return Vector3(x+b.x, y+b.y, z+b.z); }
        Vector3 operator-(const Vector3 &b) const {return Vector3(x-b.x, y-b.y, z-b.z); }
        Vector3 operator*(double b) const {return Vector3(x*b, y*b, z*b); }
        Vector3 mul(const Vector3 &b) const {return Vector3( x*b.x, y*b.y, z*b.z); }
        Vector3 cross(const Vector3 &b) const{return Vector3( y*b.z- z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }
        double dot(const Vector3 &b) const {return x*b.x+ y*b.y+ z*b.z; }
        Vector3& norm(){return *this = *this * (1/sqrt( x* x+ y* y+ z* z));}
};

struct Ray{
    Vector3 origin, dir;
    Ray(Vector3 origin_, Vector3 dir_): origin(origin_), dir(dir_) {} 
    };

struct Sphere{
    double rad; //radius
    Vector3 pos, emis, color; //pos is the center of the sphere
    Refl_t refl;
    Sphere(double rad_, Vector3 pos_, Vector3 emis_, Vector3 color_, Refl_t refl_):
        rad(rad_), pos(pos_), emis(emis_), color(color_), refl(refl_) {}
    double intersect(const Ray &r) const{
        Vector3 op = pos - r.origin;
        //solve the quadratic equation
        //P(t) = O + tD is the point along the ray
        //Equation of sphere: (P-C).(P-C) -r^2 = 0 
        //so we can get (O+tD - C). (O+tD - C) - r^2 = 0
        // => t^2(D.D) - 2tD(O-C) + (O-C).(O-C)
        // a = D.D, here a=1 because direction is unit vector, we've normalized it
        // b = 2 * D(O-C)
        // C = (O-C).(O-C) - r^2
        double t;
        double b = op.dot(r.dir); // we use 1/2 b from the above equation because fewer compuation for following det
        double det = b*b - op.dot(op) + rad*rad;// det = (B^2-4ac)/4a^2
        //since a = 1 and B = 2b => ((2b)^2-4ac)/4
        // => 4b^2 -4c /4 => b^2 - c => b^2 - (op).(op) + rad^2
        
        if(det < 0) return 0; //ray misses the sphere
        else det = sqrt(det);
        double eps = 1e-4;
        return (t=b-det) > eps ? t : ((t=b+det) > eps ? t : 0); // get the smallest positive t, the closer hit point to the ray origin
        
    }
};
//define scene: radius, position, emission, color, material
Sphere spheres[] = {
    Sphere(1e5, Vector3(1e5+1, 40.8, 81.6),   Vector3(), Vector3(.75,.25,.25), DIFF), //Left
    Sphere(1e5, Vector3(-1e5+99, 40.8, 81.6), Vector3(), Vector3(.25,.25,.75), DIFF), //Right
    Sphere(1e5, Vector3(50, 40.8, 1e5),       Vector3(), Vector3(.75,.75,.75), DIFF), //Back
    Sphere(1e5, Vector3(50, 40.8, -1e5+170),  Vector3(), Vector3(),            DIFF), //Front
    Sphere(1e5, Vector3(50, 1e5, 81.6),       Vector3(), Vector3(.75,.75,.75), DIFF), //Bottom
    Sphere(1e5, Vector3(50, -1e5+81.6, 81.6), Vector3(), Vector3(.75,.75,.75), DIFF), //Top
    Sphere(16.5, Vector3(27, 16.5, 47),       Vector3(), Vector3(1,1,1)*.999,  SPEC), //Mirror
    Sphere(16.5, Vector3(73, 16.5, 78),       Vector3(), Vector3(1,1,1)*.999,  REFR), //Glass
    Sphere(1.5, Vector3(50, 81.6-16.5, 81.6), Vector3(4,4,4)*100, Vector3(),   DIFF), //Light
};

int numSpheres = sizeof(spheres) / sizeof(Sphere);

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1: x; }
inline int toInt(double x){ return int(pow(clamp(x), 1/2.2)*255+.5); }
inline bool intersect(const Ray &r, double &t, int &id){
    double n = sizeof(spheres) / sizeof(Sphere);
    double d;
    double inf = 1e20;
    t = 1e20;
    for(int i=int(n);i--;){
        if((d = spheres[i].intersect(r))>0 && d<t){
            t = d;
            id = i;
        }
    }
    return t < inf;
}

Vector3 radiance(const Ray &r, int depth, unsigned short * Xi){
    double t;
    int id=0;
    if(!intersect(r, t, id)) return Vector3(); //if the ray miss, return black as color
    const Sphere &obj = spheres[id];

    if(depth>10) return Vector3();

    Vector3 intersect_pt = r.origin + r.dir*t;
    Vector3 n = (intersect_pt - obj.pos).norm();
    Vector3 normal = n.dot(r.dir) < 0 ? n : n*-1;
    Vector3 obj_color = obj.color;

}

int main(int argc, char *argv[]){
    //Setup image
    int width=512, height=384;
    int samples = argc==2 ? atoi(argv[1])/4 : 1;
    Vector3 r;
    Vector3 *image = new Vector3[width*height];
    //Setup camera
    //camera position and direction
    Ray cam(Vector3(50,52,295.6), Vector3(0, -0.42612, -1).norm());
    //0.5135 is the field of view angle
    Vector3 cx = Vector3(width*.5135/height, 0, 0);
    Vector3 cy = cx.cross(cam.dir).norm()*.5135;

    #pragma omp parallel for
    for(int y=0;y<height;y++){
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples*4, 100.*y/(height-1));
        unsigned short Xi[3] = {0, 0, y*y*y};
        for(unsigned int x=0;x<width;x++){
            int idx = (height-y-1)*width+x;  //calculate index for the pixel
            //for each pixel we do 2x2 subsamples
            for(int sy=0;sy<2;sy++){
                for(int sx=0;sx<2;sx++){
                    for(int s=0;s<samples;s++){
                        double r1 = 2*erand48(1);
                        double r2 = 2*erand48(1);
                        double dx = r1 < 1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
                        double dy = r2 < 1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
                        Vector3 d = cx*(((sx + 0.5 + dx)/2 + x)/width - 0.5) + 
                                    cy*(((sy + 0.5 + dy)/2 + y)/height - 0.5) + cam.dir;

                        r = r + radiance(Ray(cam.origin + d*140, d.norm()), 0, Xi)*(1./samples);
                    }
                    image[idx] = image[idx] + Vector3(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                    
                }
            }
        }

        FILE *f = fopen("image.ppm", "w");
        fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
        for(int i=0; i<width*height;i++){
            fprintf(f, "%d %d %d ", toInt(image[i].x), toInt(image[i].y), toInt(image[i].z));
        }


    }

    std::cout << "hello world!" << std::endl;
    return 0;
}
    

