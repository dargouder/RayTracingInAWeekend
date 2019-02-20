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

class XZRect : public Hitable {
public:
    XZRect() {

    }

    XZRect(float _x0, float _x1, float _z0, float _z1, float _k, std::unique_ptr<Material>
    mat) :
            x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(std::move(mat)) {

    }

    virtual bool hit(const Ray &r, float t0, float t1, HitRecord &rec) const;

    std::unique_ptr<Material> mp;
    float x0, x1, z0, z1, k;
};

bool XZRect::hit(const Ray &r, float t0, float t1, HitRecord &rec) const {
    float t = (k - r.origin().y()) / r.direction().y();
    if (t < t0 || t > t1) {
        return false;
    }

    float x = r.origin().x() + t * r.direction().x();
    float z = r.origin().z() + t * r.direction().z();

    if (x < x0 || x > x1 || z < z0 || z > z1) {
        return false;
    }

    rec.t = t;
    rec.mat_ptr = mp.get();
    rec.p = r.point_at_parameter(t);
    rec.normal = Vec3(0, 1, 0);
    return true;
}
class YZRect : public Hitable {
public:
    YZRect() {

    }

    YZRect(float _y0, float _y1, float _z0, float _z1, float _k, Material *mat) :
            y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {}

    virtual bool hit(const Ray &r, float t0, float t1, HitRecord &rec) const;

    std::unique_ptr<Material> mp;
    float y0, y1, z0, z1, k;
};


bool YZRect::hit(const Ray &r, float t0, float t1, HitRecord &rec) const {
    float t = (k - r.origin().y()) / r.direction().y();
    if (t < t0 || t > t1) {
        return false;
    }

    float y = r.origin().y() + t * r.direction().y();
    float z = r.origin().z() + t * r.direction().z();

    if (y < y0 || y > y1 || z < z0 || z > z1) {
        return false;
    }

    rec.t = t;
    rec.mat_ptr = mp.get();
    rec.p = r.point_at_parameter(t);
    rec.normal = Vec3(1, 0, 0);
    return true;
}

class FlipNormals : public Hitable
{
public:

    FlipNormals(Hitable *p) :
    shape(p)
    {}

    virtual bool hit(const Ray &r, float t_min, float t_max, HitRecord &rec) const {
        if(shape->hit(r, t_min, t_max, rec))
        {
            rec.normal = -rec.normal;
            return true;
        }
        else
        {
            return false;
        }
    }

    Hitable *shape;
};

