#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

    // Part 2, Task 4:
    // compute position and direction of ray from the input sensor sample coordinate.
    // Note: use rndR and rndTheta to uniformly sample a unit disk.
    
    Vector3D pLen(lensRadius * sqrt(rndR) * cos(rndTheta), lensRadius * sqrt(rndR) * sin(rndTheta), 0);
    float u = 2 * tan(hFov * PI / 180 / 2) * (x - 0.5);
    float v = 2 * tan(vFov * PI / 180 / 2) * (y - 0.5);

    Vector3D direction(u, v, -1);
    direction = direction - pLen / focalDistance;
    Ray r(pos + c2w * pLen, (c2w * direction).unit());
    r.min_t = fabs(nClip / direction.unit().z);
    r.max_t = fabs(fClip / direction.unit().z);
    return r;
}


} // namespace CGL
