// ShaderToy Ray Tracer for a Random Spheres Scene
#ifdef GL_ES
precision mediump float;
#endif

// iResolution, iTime, etc. are provided by ShaderToy

#define MAX_DEPTH 30
#define PI 3.14159265359

// Simple pseudo-random generator based on a 2D seed
float rand(in vec2 co) {
    return fract(sin(dot(co, vec2(12.9898,78.233))) * 43758.5453);
}

// Returns a random point in the unit sphere (using spherical coords)
vec3 random_in_unit_sphere(vec2 seed) {
    float u = rand(seed);
    float v = rand(seed+vec2(1.0,1.0));
    float theta = 2.0 * PI * u;
    float phi = acos(2.0 * v - 1.0);
    float r = pow(rand(seed+vec2(2.0,2.0)), 1.0/3.0);
    return r * vec3(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));
}



bool refracte(in vec3 v, in vec3 n, float ni_over_nt, out vec3 refracted) {
    float dt = dot(v, n);
    float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0 - dt * dt);
    if(discriminant > 0.0) {
        refracted = ni_over_nt * (v - n * dt) - n * sqrt(discriminant);
        return true;
    }
    return false;
}

float schlick(float cosine, float ref_idx) {
    float r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0*r0;
    return r0 + (1.0 - r0) * pow((1.0 - cosine), 5.0);
}

// ---------------------
// Ray & Sphere Definitions
// ---------------------
struct Ray {
    vec3 orig;
    vec3 dir;
};

vec3 ray_at(Ray r, float t) {
    return r.orig + t * r.dir;
}

// Material types: 0 = diffuse, 1 = metal, 2 = dielectric
struct Sphere {
    vec3 center;
    float radius;
    int matType;
    vec3 albedo;
   
    float ref_idx;
};

struct HitRecord {
    float t;
    vec3 p;
    vec3 normal;
    int matType;
    vec3 albedo;
    float fuzz;
    float ref_idx;
};

// Standard sphere intersection
bool hit_sphere(Sphere sph, Ray r, float t_min, float t_max, out HitRecord rec) {
    vec3 oc = r.orig - sph.center;
    float a = dot(r.dir, r.dir);
    float half_b = dot(oc, r.dir);
    float c = dot(oc, oc) - sph.radius * sph.radius;
    float discriminant = half_b * half_b - a * c;
    if(discriminant > 0.0) {
        float sqrtd = sqrt(discriminant);
        float root = (-half_b - sqrtd) / a;
        if(root < t_max && root > t_min) {
            rec.t = root;
            rec.p = ray_at(r, rec.t);
            rec.normal = (rec.p - sph.center) / sph.radius;
            rec.matType = sph.matType;
            rec.albedo = sph.albedo;
            
            rec.ref_idx = sph.ref_idx;
            return true;
        }
        root = (-half_b + sqrtd) / a;
        if(root < t_max && root > t_min) {
            rec.t = root;
            rec.p = ray_at(r, rec.t);
            rec.normal = (rec.p - sph.center) / sph.radius;
            rec.matType = sph.matType;
            rec.albedo = sph.albedo;
           
            rec.ref_idx = sph.ref_idx;
            return true;
        }
    }
    return false;
}

// ---------------------
// Scene Definition
// ---------------------
bool hit_world(Ray r, float t_min, float t_max, out HitRecord rec) {
    HitRecord tempRec;
    bool hitAnything = false;
    float closestSoFar = t_max;
    
    // Ground sphere
    Sphere ground;
    ground.center = vec3(0.0, -1000.0, 0.0);
    ground.radius = 1000.0;
    ground.matType = 0;
    ground.albedo = vec3(0.5);
    
    ground.ref_idx = 1.0;
    if(hit_sphere(ground, r, t_min, closestSoFar, tempRec)) {
        hitAnything = true;
        closestSoFar = tempRec.t;
        rec = tempRec;
    }
    
    // Three main spheres
    Sphere s1; // dielectric
    s1.center = vec3(0.0, 1.0, 0.0);
    s1.radius = 1.0;
    s1.matType = 2;
    s1.albedo = vec3(0.8,0.8,0.8);
    
    s1.ref_idx = 1.5;
    if(hit_sphere(s1, r, t_min, closestSoFar, tempRec)) {
        hitAnything = true;
        closestSoFar = tempRec.t;
        rec = tempRec;
    }
    
    Sphere s2; // diffuse
    s2.center = vec3(-4.0, 1.0, 0.0);
    s2.radius = 1.0;
    s2.matType = 0;
    s2.albedo = vec3(0.4, 0.2, 0.1);
    
    s2.ref_idx = 1.0;
    if(hit_sphere(s2, r, t_min, closestSoFar, tempRec)) {
        hitAnything = true;
        closestSoFar = tempRec.t;
        rec = tempRec;
    }
    
    Sphere s3; // metal
    s3.center = vec3(4.0, 1.0, 0.0);
    s3.radius = 1.0;
    s3.matType = 1;
    s3.albedo = vec3(0.7, 0.6, 0.5);

    s3.ref_idx = 1.0;
    if(hit_sphere(s3, r, t_min, closestSoFar, tempRec)) {
        hitAnything = true;
        closestSoFar = tempRec.t;
        rec = tempRec;
    }
    
    // Small spheres in a grid from a=-11 to 10 and b=-11 to 10
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            vec2 seed = vec2(float(a), float(b));
            float choose_mat = rand(seed);
            vec3 center = vec3(float(a) + 0.9 * rand(seed+vec2(1.0,2.0)),
                               0.2,
                               float(b) + 0.9 * rand(seed+vec2(3.0,4.0)));
            // Avoid overlap with big sphere at (4,0.2,0)
            if(length(center - vec3(4.0, 0.2, 0.0)) > 0.9) {
                Sphere s;
                s.center = center;
                s.radius = 0.2;
                if (choose_mat < 0.8) {
                    // Diffuse
                    s.matType = 0;
                    s.albedo = vec3(rand(seed+vec2(5.0,6.0))*rand(seed+vec2(7.0,8.0)),
                                    rand(seed+vec2(9.0,10.0))*rand(seed+vec2(11.0,12.0)),
                                    rand(seed+vec2(13.0,14.0))*rand(seed+vec2(15.0,16.0)));
                   
                    s.ref_idx = 1.0;
                } else if (choose_mat < 0.95) {
                    // Metal
                    s.matType = 1;
                    s.albedo = vec3(rand(seed+vec2(17.0,18.0)) * 0.5 + 0.5,
                                    rand(seed+vec2(19.0,20.0)) * 0.5 + 0.5,
                                    rand(seed+vec2(21.0,22.0)) * 0.5 + 0.5);
                  
                    s.ref_idx = 1.0;
                } else {
                    // Glass (dielectric)
                    s.matType = 2;
                    s.albedo = vec3(1.0);
                 
                    s.ref_idx = 1.5;
                }
                if(hit_sphere(s, r, t_min, closestSoFar, tempRec)) {
                    hitAnything = true;
                    closestSoFar = tempRec.t;
                    rec = tempRec;
                }
            }
        }
    }
    
    return hitAnything;
}

// ---------------------
// Ray Color: Traces the ray and computes attenuation
// ---------------------
vec3 ray_color(Ray r, vec2 seed) {
    vec3 col = vec3(0.0);
    vec3 attenuation = vec3(1.0);
    for (int i = 0; i < MAX_DEPTH; i++) {
        HitRecord rec;
        if (hit_world(r, 0.001, 10000.0, rec)) {
            Ray scattered;
            if (rec.matType == 0) {
                // Diffuse: scatter in random direction in hemisphere
                vec3 scatterDir = rec.normal ;
                if(length(scatterDir) < 0.001)
                    scatterDir = rec.normal;
                scattered.orig = rec.p;
                scattered.dir = scatterDir;
                attenuation *= rec.albedo;
                r = scattered;
            } else if (rec.matType == 1) {
                // Metal: reflect with fuzz
                vec3 reflected = reflect(normalize(r.dir), rec.normal);
                scattered.orig = rec.p;
                scattered.dir = reflected + rec.fuzz ;
                if(dot(scattered.dir, rec.normal) <= 0.0) {
                    attenuation = vec3(0.0);
                    break;
                }
                attenuation *= rec.albedo;
                r = scattered;
            } else if (rec.matType == 2) {
                // Dielectric (glass)
                attenuation *= vec3(1.0);
                float ref_ratio = (dot(r.dir, rec.normal) < 0.0) ? (1.0 / rec.ref_idx) : rec.ref_idx;
                vec3 unitDir = normalize(r.dir);
                float cosTheta = min(dot(-unitDir, rec.normal), 1.0);
                float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
                vec3 direction;
                if(ref_ratio * sinTheta > 1.0 || schlick(cosTheta, ref_ratio) > rand(seed+vec2(float(i)))) {
                    direction = reflect(unitDir, rec.normal);
                } else {
                    vec3 refracted;
                    refracte(unitDir, rec.normal, ref_ratio, refracted);
                    direction = refracted;
                }
                scattered.orig = rec.p;
                scattered.dir = direction;
                r = scattered;
            }
        } else {
            // Background gradient
            vec3 unitDir = normalize(r.dir);
            float t = 0.5 * (unitDir.y + 1.0);
            vec3 background = mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            col = attenuation * background;
            break;
        }
    }
    return col;
}

// ---------------------
// Camera Setup
// ---------------------
struct Camera {
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    float lens_radius;
};

Camera camera_setup(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov, float aspect, float focus_dist) {
    Camera cam;
    float theta = radians(vfov);
    float h = tan(theta/2.0);
    float viewport_height = 2.0 * h;
    float viewport_width = aspect * viewport_height;
    
    cam.w = normalize(lookfrom - lookat);
    cam.u = normalize(cross(vup, cam.w));
    cam.v = cross(cam.w, cam.u);
    
    cam.origin = lookfrom;
    cam.horizontal = focus_dist * viewport_width * cam.u;
    cam.vertical = focus_dist * viewport_height * cam.v;
    cam.lower_left_corner = cam.origin - cam.horizontal * 0.5 - cam.vertical * 0.5 - focus_dist * cam.w;
 
    return cam;
}

Ray get_ray(Camera cam, vec2 uv, vec2 seed) {
    // Depth-of-field: random point in unit disk
    float r = rand(seed);
    float angle = 2.0 * PI * rand(seed+vec2(1.0,1.0));
    vec2 rd = cam.lens_radius * r * vec2(cos(angle), sin(angle));
    vec3 offset = cam.u * rd.x + cam.v * rd.y;
    Ray rayo;
    rayo.orig = cam.origin + offset;
    rayo.dir = normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset);
    return rayo;
}

// ---------------------
// Main Shader Function
// ---------------------
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalize pixel coordinates
    vec2 uv = fragCoord / iResolution.xy;
    
    // Setup camera (matching your C++ code)
    vec3 lookfrom = vec3(13.0, 2.0, 3.0);
    vec3 lookat   = vec3(0.0, 0.0, 0.0);
    vec3 vup      = vec3(0.0, 1.0, 0.0);
    float vfov = 20.0;
    float aspect = iResolution.x / iResolution.y;
    float aperture = 0.6;
    float focus_dist = 100.0;
    Camera cam = camera_setup(lookfrom, lookat, vup, vfov, aspect, focus_dist);
    
    vec3 col = vec3(0.0);
    int samples = 3; // Increase for better quality if performance allows
    for (int s = 0; s < samples; s++) {
        // Anti-alias by jittering pixel coordinates
        vec2 jitter = vec2(rand(fragCoord + vec2(float(s))), rand(fragCoord + vec2(float(s)+1.0)));
        vec2 uvj = (fragCoord + jitter) / iResolution.xy;
        Ray r = get_ray(cam, uvj, fragCoord + vec2(float(s)));
        col += ray_color(r, fragCoord + vec2(float(s)));
    }
    col /= float(samples);
    
    // Gamma correction
    col = sqrt(col);
    
    fragColor = vec4(col, 1.0);}
