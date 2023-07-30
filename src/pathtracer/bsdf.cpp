#include "bsdf.h"
#include "bsdf.h"
#include "bsdf.h"
#include "bsdf.h"

#include "application/visual_debugger.h"

#include <algorithm>
#include <iostream>
#include <utility>


using std::max;
using std::min;
using std::swap;

namespace CGL {

/**
 * This function creates a object space (basis vectors) from the normal vector
 */
void make_coord_space(Matrix3x3 &o2w, const Vector3D n) {

  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z))
    h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z))
    h.y = 1.0;
  else
    h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

/**
 * Evaluate diffuse lambertian BSDF.
 * Given incident light direction wi and outgoing light direction wo. Note
 * that both wi and wo are defined in the local coordinate system at the
 * point of intersection.
 * \param wo outgoing light direction in local space of point of intersection
 * \param wi incident light direction in local space of point of intersection
 * \return reflectance in the given incident/outgoing directions
 */
Vector3D DiffuseBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO (Part 3.1):
  // This function takes in both wo and wi and returns the evaluation of
  // the BSDF for those two directions.
    return reflectance / (2*PI);
}

/**
 * Evalutate diffuse lambertian BSDF.
 */
Vector3D DiffuseBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf) {
  // TODO (Part 3.1):
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).
  // You can use the `f` function. The reference solution only takes two lines.

    *wi = sampler.get_sample(pdf);
    return f(wo, *wi); 
}

void DiffuseBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Diffuse BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}


Vector3D FogBSDF::f(const Vector3D wo, const Vector3D wi) {
    return reflectance/4/PI;
}

Vector3D FogBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
    //*wi = sampler.get_sample();
    //double cos_wi_rd = dot(wi->unit(), wo.unit());
    //*pdf = 1 / (4 * PI) * (1 - phase_factor* phase_factor) / (pow(1 + phase_factor * phase_factor - 2 * phase_factor * cos_wi_rd, 3 / 2));
    //return f(wo, *wi);
    double Xi1 = random_uniform();
    double Xi2 = random_uniform();

    double r = sqrt(Xi1);
    double theta = 2. * PI * Xi2;
    if (coin_flip(0.5)) *wi = Vector3D(sqrt(1 - Xi1), r * cos(theta), r * sin(theta));
    else *wi = Vector3D(-sqrt(1 - Xi1), r * cos(theta), r * sin(theta));
    *pdf = sqrt(1 - Xi1) / PI/2;
    return f(wo, *wi);
}

double FogBSDF::get_HG(const Vector3D wo, Vector3D wi)
{
    double cos_wi_rd = dot(wi.unit(), wo.unit());
    return 1 / (4 * PI) * (1 - phase_factor * phase_factor) / (pow(1 + phase_factor * phase_factor - 2 * phase_factor * cos_wi_rd, 3 / 2)); 
}



/**
 * Evalutate Emission BSDF (Light Source)
 */
Vector3D EmissionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

/**
 * Evalutate Emission BSDF (Light Source)
 */
Vector3D EmissionBSDF::sample_f(const Vector3D wo, Vector3D *wi, double *pdf) {
  *pdf = 1.0 / PI;
  *wi = sampler.get_sample(pdf);
  return Vector3D();
}

void EmissionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Emission BSDF"))
  {
    DragDouble3("Radiance", &radiance[0], 0.005);
    ImGui::TreePop();
  }
}

} // namespace CGL
