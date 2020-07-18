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
    bool isLight() const{
        return emis.x<=0 && emis.y<=0 && emis.z<=0;
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

Vector3 radiance(const Ray &r, int depth, unsigned short * Xi, int E=1){
    double t;
    int id=0;
    if(!intersect(r, t, id)) return Vector3(); //if the ray miss, return black as color
    const Sphere &obj = spheres[id];

    if(depth>10) return Vector3();

    Vector3 intersect_pt = r.origin + r.dir*t;
    Vector3 n = (intersect_pt - obj.pos).norm();
    Vector3 normal = n.dot(r.dir) < 0 ? n : n*-1;
    Vector3 f = obj.color;

    //russian roulette
    double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z;
    if(++depth>5||!p)
    {
        if(erand48(1)<p)
            f = f*(1/p);
        else 
            return obj.emis*E;
    }

    if(obj.refl == DIFF){  //Ideal DIFFUSE reflection
        double r1 = 2*M_PI*erand48(1); //get random angle
        double r2 = erand48(1), r2s=sqrt(r2); 
        Vector3 w = normal;
        Vector3 u = ((fabs(w.x) > .1? Vector3(0,1):Vector3(1)).cross(w)).norm();//u is perpendicular to w
        Vector3 v = w.cross(u); //v is perpendicular to w and u
        Vector3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();

        //Loop over all lights(explicit lighting)
        Vector3 e;
        for(int i=0;i<numSpheres;i++){
            const Sphere &s = spheres[i];
            if(!s.isLight()) continue; //skip the sphere that is not light

            //Create random direction 
            Vector3 sw = s.pos - intersect_pt;
            Vector3 su = ((fabs(sw.x) > .1? Vector3(0,1):Vector3(1)).cross(sw)).norm();
            Vector3 sv = sw.cross(su);
            //TODO: do the math again
            double cos_a_max = sqrt(1-s.rad*s.rad/(intersect_pt-s.pos).dot(intersect_pt-s.pos));
            double eps1 = erand48(1), eps2 =erand48(1);
            double cos_a = 1-eps1+eps1*cos_a_max;
            double sin_a = sqrt(1-cos_a*cos_a);
            double phi = 2*M_PI*eps2;
            Vector3 l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
            l.norm();

            //shoot shadow ray
            if(intersect(Ray(intersect_pt,l), t, id) && id==i){
                double omega = 2*M_PI*(1-cos_a_max);
                e = e +f.mul(s.emis*l.dot(normal)*omega)*M_1_PI;
            }
        }
        return obj.emis*E + e + f.mul(radiance(Ray(intersect_pt, d), depth, Xi, 0));
    }
    else if(obj.refl == SPEC){
        return obj.emis + f.mul(radiance(Ray(intersect_pt,r.dir-n*2*n.dot(r.dir)), depth, Xi));
    }

    Ray reflRay(intersect_pt, r.dir-n*2*n.dot(r.dir));
    bool into = n.dot(normal) > 0;
    double nc = 1, nt = 1.5;
    double nnt = into?nc/nt:nt/nc;
    double ddn = r.dir.dot(normal), cos2t; 

    if((cos2t = 1-nnt*nnt*(1-ddn*ddn))<0)
        return obj.emis + f.mul(radiance(reflRay, depth, Xi));
    Vector3 tdir = (r.dir * nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c= 1-(into?-ddn:tdir.dot(n));
    double Re = R0+(1-R0)*c*c*c*c*c, Tr=1-Re, P=.25+.5*Re, RP=Re/P, TP=Tr/(1-P);
    return obj.emis + f.mul(depth>2 ? (erand48(1)<P ?
    radiance(reflRay, depth, Xi)*RP: radiance(Ray(intersect_pt, tdir), depth, Xi)*TP):
    radiance(reflRay, depth, Xi)* Re+radiance(Ray(intersect_pt, tdir), depth, Xi)*Tr);

}


int main(int argc, char *argv[]){
    //const double M_PI_ = 3.1415926535;
    //const double M_1_PI = 1.0/M_PI;
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
    

