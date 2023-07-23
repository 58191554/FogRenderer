#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

    // TODO (Part 2.2):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bouding box within the range given by
    // t0, t1, update t0 and t1 with the new intersection times.
    
    // yz-planes

    double tx0 = (min[0] - r.o[0]) / r.d[0];
    double tx1 = (max[0] - r.o[0]) / r.d[0];
    if (tx1 < tx0) std::swap(tx0, tx1);

    // xz-planes

    double ty0 = (min[1] - r.o[1]) / r.d[1];
    double ty1 = (max[1] - r.o[1]) / r.d[1];
    if (ty1 < ty0) std::swap(ty0, ty1);

    // xy-planes
    double tz0 = (min[2] - r.o[2]) / r.d[2];
    double tz1 = (max[2] - r.o[2]) / r.d[2];
    if (tz1 < tz0) std::swap(tz0, tz1);

    // Check rectangle intersection along a single z-slice first 
    double t_min = std::max(tx0, ty0);
    double t_max = std::min(tx1, ty1);
    if (t_min > t_max) return false;

    // If still a contender, check intersection points along z planes
    if (tz0 > t_max || tz1 < t_min) return false;
    t_min = std::max(tz0, t_min);
    t_max = std::min(tz1, t_max);

    // Check that entrance (t_min) and exit (t_max) times make sense.
    if (t_min > t_max) return false;

    // Check that there is overlap between the ray interval and the box interesection interval.
    if (!((t1 - t_min) >= 0 && (t_max - t0) >= 0)) return false;

    // All cases checked so return that here was indeed an intersection with a box.
    return true;
}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
