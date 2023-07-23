#include "pathtracer.h"
#include "pathtracer.h"
#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
    // Estimate the lighting from this intersection coming directly from a light.
    // For this function, sample uniformly in a hemisphere.

    // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
    // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

    // make a coordinate system for a hit point
    // with N aligned with the Z direction.
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();

    // w_out points towards the source of the ray (e.g.,
    // toward the camera if this is a primary ray)
    const Vector3D hit_p = r.o + r.d * isect.t;
    const Vector3D w_out = w2o * (-r.d);

    // This is the same number of total samples as
    // estimate_direct_lighting_importance (outside of delta lights). We keep the
    // same number of samples for clarity of comparison.
    int num_samples = scene->lights.size() * ns_area_light;
    Vector3D L_out(0.0, 0.0, 0.0);

    // TODO (Part 3): Write your sampling loop here
    // TODO BEFORE YOU BEGIN
    // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 

    for (int i = 0; i < num_samples; i++) {
        Vector3D w_in_i;
        double pdf;
        Vector3D ref_rad = isect.bsdf->sample_f(w_out, &w_in_i, &pdf);
        float cos_theta = dot(w_in_i, Vector3D(0, 0, 1));
        Ray r_in_i = Ray(hit_p, o2w * w_in_i, 1);
        r_in_i.min_t = EPS_F;
        Intersection isect_i;
        if (cos_theta > 0 && bvh->intersect(r_in_i, &isect_i)) {
            Vector3D L_i = isect_i.bsdf->get_emission();
            L_out += L_i * 2*PI*ref_rad * cos_theta/pdf;
        }
        L_out /= num_samples;
    }
    return L_out;
}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
    // Estimate the lighting from this intersection coming directly from a light.
    // To implement importance sampling, sample only from lights, not uniformly in
    // a hemisphere.

    // make a coordinate system for a hit point
    // with N aligned with the Z direction.
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();

    // w_out points towards the source of the ray (e.g.,
    // toward the camera if this is a primary ray)
    const Vector3D hit_p = r.o + r.d * isect.t;
    const Vector3D w_out = w2o * (-r.d);
    Vector3D L_out, w_in;
    double dist, pdf;
    int num_sample;
    for (auto light = scene->lights.begin(); light != scene->lights.end(); light++) {
        if ((*light)->is_delta_light()) num_sample = 1;
        else num_sample = ns_area_light;
        for (int i = 0; i < num_sample; ++i) {
            Vector3D L_i = (*light)->sample_L(hit_p, &w_in, &dist, &pdf);
            Ray r_i(hit_p, w_in.unit());
            r_i.min_t = EPS_F; r_i.max_t = dist;
            Intersection isect_i;
            w_in = w2o * w_in;
            if (w_in[2] < 0) continue;
            if (!bvh->intersect(r_i, &isect_i)) {
                Vector3D dL = isect.bsdf->f(w_out, w_in) * L_i * w_in[2]*PI / pdf;
                L_out += (1.0 / num_sample) * dL;
            }
        }
    }
    return L_out;
}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
    // TODO: Part 3, Task 2
    // Returns the light that results from no bounces of light
    return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
    if (direct_hemisphere_sample) {
        return estimate_direct_lighting_hemisphere(r, isect);
    }
    else {
        return estimate_direct_lighting_importance(r, isect);
    }
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();

    Vector3D hit_p = r.o + r.d * isect.t;
    Vector3D w_out = w2o * (-r.d);

    Vector3D L_out;

    // TODO: Part 4, Task 2
    // Returns the one bounce radiance + radiance from extra bounces at this point.
    // Should be called recursively to simulate extra bounces.
    Vector3D w_in;
    double pdf;
    Vector3D f = isect.bsdf->sample_f(w_out, &w_in, &pdf);
    Ray new_r(hit_p, o2w * w_in);
    new_r.depth = r.depth - 1;
    new_r.min_t = EPS_D;
    new_r.max_t = INF_D;
    Intersection new_isect;
    bool flag = bvh->intersect(new_r, &new_isect);
    double random_prob = 0.35;
    if (flag && r.depth > 0 && coin_flip(random_prob))return one_bounce_radiance(r, isect) + f * at_least_one_bounce_radiance(new_r, new_isect) * w_in[2] / pdf / random_prob;
    else return one_bounce_radiance(r, isect);
}

Vector3D PathTracer::fog_bounce_radiance(const Ray& r, const SceneObjects::Intersection& isect)
{
    // TODO: importance sampling uniformly sample (isotropic)
    return estimate_fog_lighting_sphere(r, isect);

}

double calculateCosine(const Vector3D& v1, const Vector3D& v2) {
    double dotProduct = dot(v1, v2);
    double v1Length = v1.norm(); 
    double v2Length = v2.norm(); 

    double cosine = dotProduct / (v1Length * v2Length); 
    return cosine;
}

Vector3D PathTracer::estimate_fog_lighting_sphere(const Ray& r, const SceneObjects::Intersection& isect)
{
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();
    double g = 0;

    // w_out points towards the source of the ray (e.g.,
    // toward the camera if this is a primary ray)
    const Vector3D hit_p = r.o + r.d * isect.t;
    const Vector3D w_out = w2o * (-r.d);
    Vector3D L_out, w_in;
    double dist, pdf;
    int num_sample;
    for (auto light = scene->lights.begin(); light != scene->lights.end(); light++) {
        if ((*light)->is_delta_light()) num_sample = 1;
        else num_sample = ns_area_light;
        for (int i = 0; i < num_sample; ++i) {
            Vector3D L_i = (*light)->sample_L(hit_p, &w_in, &dist, &pdf);
            Ray r_i(hit_p, w_in.unit());
            r_i.min_t = EPS_F; r_i.max_t = dist;
            Intersection isect_i;
            w_in = w2o * w_in;
            if (w_in[2] < 0) continue;
            if (!bvh->intersect(r_i, &isect_i)) {
                double HG = 1 / (4 * PI) * (1 - g * g) / (pow(1 + g * g - 2 * g * calculateCosine(r.d, w_in), 3 / 2));
                Vector3D dL = isect.bsdf->f(w_out, w_in) / 2 * L_i * w_in[2] / HG / pdf;
                //Vector3D dL = isect.bsdf->f(w_out, w_in)/2 * L_i * w_in[2] * PI / pdf;
                L_out += (1.0 / num_sample) * dL;
            }
        }
    }
    return L_out;
}

double sampleInverseExponential(double lambda) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    double u = dis(gen);
    return -log(1.0 - u) / lambda;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
    Intersection isect;
    Vector3D L_out;
    UniformSphereSampler3D sphere_sampler;

    // You will extend this in assignment 3-2.

    // If no intersection occurs, we simply return black.
    // This changes if you implement hemispherical lighting for extra credit.

    bool FOG_FLAG = true;
    if (FOG_FLAG) {
        Intersection fog_isect;
        double LAMBDA = 0.02;
        // randomly intersect fog partical 
        fog_isect.t = sampleInverseExponential(LAMBDA) + r.min_t;
        //fog_isect.t = random_uniform() * (r.max_t - r.min_t) + r.min_t;
        bvh->intersect(r, &isect);
        if (!bvh->intersect(r, &isect) && r.max_t < fog_isect.t)
            return envLight ? envLight->sample_dir(r) : L_out;

        if (isect.t > fog_isect.t) {
            fog_isect.n = sphere_sampler.get_sample();
            fog_isect.bsdf = new DiffuseBSDF(Vector3D(1, 1, 1));
            L_out += estimate_fog_lighting_sphere(r, fog_isect);
            return L_out;
        }

        //L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);
        Vector3D L_r;
        L_r = zero_bounce_radiance(r, isect);
        L_r += at_least_one_bounce_radiance(r, isect);
        double sigma = 0.3;
        L_out = L_r * exp(-sigma * isect.t);
        return L_out;
    }
    else {
        if (!bvh->intersect(r, &isect))
            return envLight ? envLight->sample_dir(r) : L_out;
        L_out = zero_bounce_radiance(r, isect);
        L_out += at_least_one_bounce_radiance(r, isect);
        return L_out;
    }
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel

  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  Vector3D radiance(0, 0, 0);
  double s1, s2 = 0;
  double illu;
  int i;
  for (i = 0; i < num_samples; i++) {
      Vector2D random_sample = gridSampler->get_sample();
      Vector2D pixel_sample = origin + random_sample;
      Ray r = camera->generate_ray(pixel_sample.x / sampleBuffer.w, pixel_sample.y / sampleBuffer.h);
      r.depth = max_ray_depth;
      Vector3D rad_i = est_radiance_global_illumination(r);
      radiance += rad_i;
      illu = rad_i.illum();
      s1 += illu;
      s2 += illu * illu;
      if (i % samplesPerBatch == 0) {
          if (1.96 * sqrt((s2 - s1 * s1 / i) / (i - 1)) / sqrt(i) < maxTolerance * (s1 / i)) break;
      }

  }
  radiance /= i;

  sampleBuffer.update_pixel(radiance, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
