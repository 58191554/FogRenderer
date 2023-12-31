#ifndef CGL_COLLADA_LIGHTINFO_H
#define CGL_COLLADA_LIGHTINFO_H

#include "CGL/vector3D.h"

#include "collada_info.h"

namespace CGL { namespace Collada {

// For some reason enum classes aren't working for me here; not sure why.
namespace LightType {
  enum T {
    NONE,
    AMBIENT,
    DIRECTIONAL,
    AREA,
    POINT,
    SPOT
  };
}

class LightInfo : public Instance {
 public:
  
  LightInfo();

  LightType::T light_type;  ///< type 

  Vector3D spectrum;        ///< spectrum 

  Vector3D position;        ///< position
  Vector3D direction;       ///< direction
  Vector3D up;              ///< up

  float falloff_deg;        ///< degree of fall off angle
  float falloff_exp;        ///< fall out exponent

  float constant_att;       ///< constant attenuation factor  
  float linear_att;         ///< linear attenuation factor
  float quadratic_att;      ///< quadratic attenuation factor
  float falloff_angle;      ///< fall of angle factor Zhen TONG added

}; // struct Light

std::ostream& operator<<( std::ostream& os, const LightInfo& light );

} // namespace Collada
} // namespace CGL

#endif // CGL_COLLADA_LIGHTINFO_H
