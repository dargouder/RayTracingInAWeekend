//
// Created by Darryl Gouder on 2019-02-18.
//
#pragma once

#include "hitable.h"

class XYRect : public Hitable {
public:
    XYRect() {

    }

    XYRect(float _x0, float _x1, float _y0, float _y1, float _k, std::unique_ptr<Material>
    mat) :
            x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(std::move(mat)) {

    }

    virtual bool hit(const Ray &r, float t0, float t1, HitRecord &rec) const;

    std::unique_ptr<Material> mp;
    float x0, x1, y0, y1, k;
};

bool XYRect::hit(const Ray &r, float t0, float t1, HitRecord &rec) const {
    float t = (k - r.origin().y()) / r.direction().y();
    if (t < t0 || t > t1) {
        return false;
    }

    float x = r.origin().x() + t * r.direction().x();
    float y = r.origin().y() + t * r.direction().y();

    if (x < x0 || x > x1 || y < y0 || y > y1) {
        return false;
    }

    rec.t = t;
    rec.mat_ptr = mp.get();
    rec.p = r.point_at_parameter(t);
    rec.normal = Vec3(0, 0, 1);
    return true;
}

//class YZRect : public Hitable {
//public:
//    YZRect() {
//
//    }
//
//    YZRect(float _y0, float _y1, float _z0, float _z1, float _k, Material *mat) :
//            y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {}
//
//    virtual bool hit(const Ray &r, float t0, float t1, HitRecord &rec) const;
//
//    Material *mp;
//    float x0, x1, z0, z1, k;
//};
//
//bool YZRect::hit(const Ray &r float t0, float t1, HitRecord &rec) const {
//
//}

