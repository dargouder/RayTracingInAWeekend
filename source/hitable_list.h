#pragma once
#include <vector>
#include <memory>
#include "hitable.h"


class HitableList : public Hitable {
public:
	HitableList() {}

	bool hit(const Ray& ray, float tmin, float tmax, HitRecord& rec) const {
		HitRecord temp_rec;

		bool hit_anything = false;

		double closest_so_far = tmax;

		for(auto& hitable_object : list) {
			if(hitable_object->hit(ray, tmin, closest_so_far,temp_rec))
			{
				hit_anything = true;
				closest_so_far = temp_rec.t;
				rec = temp_rec;
			}
		}

		return hit_anything;
	}

  Vec3 generateSampleOnSurface() const
  {
    return Vec3(0.0f, 0.0f, 0.0f);
  }

	std::vector<std::unique_ptr<Hitable>> list;
};