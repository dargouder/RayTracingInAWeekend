#include "reflection.h"

Vec3 BxDF::Sample_f(const Vec3 &wo, Vec3 &wi, float &pdf) const {
  // Cosine-sample the hemisphere, flipping the direction if necessary
  // direction is in the hemisphere around (0, 0, 1)
  wi = CosineSampleHemispherePBRT(RAND(), RAND());
  if (wo.z < 0) {
    // wo was pointing in the opposite hemisphere
    wi.z *= -1;
  }
  pdf = Pdf(wo, wi);
  return f(wo, wi);
}

float BxDF::Pdf(const Vec3 &wo, const Vec3 &wi) const {
  // if they're on the same hemisphere, the pdf is | n . wi |
  // Since our local normal is (0, 0, 1), we can use the utility BXDF functions
  return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * INV_PI : 0;
}
