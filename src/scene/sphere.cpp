#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.


  return true;

}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
    double a, b, c, t, t1, t2, discriminant;

    bool valid_intersection = false;

    a = dot(r.d, r.d);
    b = dot(2 * (r.o - o), r.d);
    c = dot((r.o - o), (r.o - o)) - r2;

    discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {

        t1 = (-b + sqrt(discriminant)) / (2 * a);
        t2 = (-b - sqrt(discriminant)) / (2 * a);

        t = t1 < t2 ? t1 : t2;

        valid_intersection = r.min_t <= t && t <= r.max_t;
    }

    else if (discriminant == 0) {
        t = -b / (2 * a);

        valid_intersection = r.min_t <= t && t <= r.max_t;
    }


    if (valid_intersection) {
        r.max_t = t;
        return true;
    }
    else {
        return false;
    }
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
    if (has_intersection(r)) {
        double a, b, c, t, t1, t2, discriminant;

        a = dot(r.d, r.d);
        b = dot(2 * (r.o - o), r.d);
        c = dot((r.o - o), (r.o - o)) - r2;

        discriminant = b * b - 4 * a * c;

        if (discriminant > 0) {
            t1 = (-b + sqrt(discriminant)) / (2 * a);
            t2 = (-b - sqrt(discriminant)) / (2 * a);
            t = t1 < t2 ? t1 : t2;
        }

        else {
            t = -b / (2 * a);
        }

        //compute surface normal;
        Vector3D intersection_point = r.o + t * r.d;
        Vector3D normal = intersection_point - o;
        normal.normalize();

        //fill in intersection data structure.
        i->t = t;
        i->primitive = this;
        i->bsdf = get_bsdf();
        i->n = normal;
        return true;
    }
    else {
        return false;
    }
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
