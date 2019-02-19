//
// Created by Darryl Gouder on 2019-02-18.
//
#pragma once
#include "hitable.h

class XZRect : public Hitable
{
public:
    XZRect()
    {

    }

    XZRect(float _x0, float _x1, float _z0, float _z1, float _k, Material *mat) :
    x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat)
    {}

    virtual bool hit(const Ray &r, float t0, float t1, HitRecord &rec) const;

    Material *mp;
    float x0, x1, z0, z1, k;
};

bool XZRect::hit(const Ray &r float t0, float t1, HitRecord &rec) const
{
float t = (k-r.origin().y()) / r.direction().y();
if(t < t0 ||)
}

class YZRect : public Hitable
{
public:
    YZRect()
    {

    }

    YZRect(float _y0, float _y1, float _z0, float _z1, float _k, Material *mat) :
    y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat)
    {}

    virtual bool hit(const Ray&r, float t0, float t1, HitRecord &rec) const;

    Material *mp;
    float x0, x1, z0, z1, k;
};

bool YZRect::hit(const Ray &r float t0, float t1, HitRecord &rec) const
{

}

