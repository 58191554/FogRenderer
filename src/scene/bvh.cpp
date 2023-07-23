#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>
#include <fstream>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size, 0);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode* BVHAccel::construct_bvh(std::vector<Primitive*>::iterator start,
    std::vector<Primitive*>::iterator end,
    size_t max_leaf_size, int cur_depth) {

    // TODO (Part 2.1):
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration. The starter code build a BVH aggregate with a
    // single leaf node (which is also the root) that encloses all the
    // primitives.


    BBox bbox;

    for (auto p = start; p != end; p++) {
        BBox bb = (*p)->get_bbox();
        bbox.expand(bb);
    }

    BVHNode* node = new BVHNode(bbox);// create the return node
    // too less primitives to be a intermediate node OR too deep to divide
    if (end - start < max_leaf_size) {
        node->start = start;
        node->end = end;
        node->l = NULL;
        node->r = NULL;
        node->depth = cur_depth + 1;
        return node;
    }

    // split into 2 subset of primitives
    Vector3D mean(0.0, 0.0, 0.0);
    for (auto p = start; p != end; p++) {
        mean += (*p)->get_bbox().centroid();
    }
    mean /= (end - start);// Get the mean cetroid position of whole set of primitives
    float min_heu_val = INFINITY;// set the initial heuristic value for comparasion;

    std::vector<Primitive*> l_ps, r_ps;
    BBox l_bb, r_bb;
    for (int i = 0; i < 3; i++) {
        std::vector<Primitive*> tmp_l_ps, tmp_r_ps;
        BBox tmp_l_bb, tmp_r_bb;
        for (auto p = start; p != end; p++) {
            if ((*p)->get_bbox().centroid()[i] <= mean[i]) {
                tmp_l_ps.push_back(*p);
                tmp_l_bb.expand((*p)->get_bbox());
            }
            else {
                tmp_r_ps.push_back(*p);
                tmp_r_bb.expand((*p)->get_bbox());
            }
        }
        // compute the heuristic of the split
        Vector3D exL = tmp_l_bb.extent;
        Vector3D exR = tmp_r_bb.extent;
        float heuristic = tmp_l_ps.size() * (exL.x * exL.y + exL.y * exL.z + exL.z * exL.x)
            + tmp_r_ps.size() * (exR.x * exR.y + exR.y * exR.z + exR.z * exR.x);
        if (heuristic < min_heu_val) {
            min_heu_val = heuristic;
            l_ps = tmp_l_ps;
            r_ps = tmp_r_ps;
            l_bb = tmp_l_bb;
            r_bb = tmp_r_bb;
        }
    }

    if (l_ps.size() == end - start || r_ps.size() == end - start) {
        if (l_ps.size() > r_ps.size()) {
            r_ps.push_back(l_ps.back());
            l_ps.pop_back();
        }
        else
        {
            l_ps.push_back(r_ps.back());
            r_ps.pop_back();
        }
    }

    //update elements in orginall array
    for (int i = 0; i < l_ps.size(); i++) {
        *(start + i) = l_ps.at(i);
    }

    for (int j = 0; j < r_ps.size(); j++) {
        *(start + l_ps.size() + j) = r_ps.at(j);
    }


    node->l = construct_bvh(start, next(start, l_ps.size()), max_leaf_size, cur_depth + 1);
    node->r = construct_bvh(next(start, l_ps.size()), end, max_leaf_size, cur_depth + 1);

    return node;
}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
    // TODO (Part 2.3):
    // Fill in the intersect function.
    // Take note that this function has a short-circuit that the
    // Intersection version cannot, since it returns as soon as it finds
    // a hit, it doesn't actually have to find the closest hit.
    if (node->bb.intersect(ray, ray.min_t, ray.max_t)) {
        if (node->isLeaf()) {
            for (auto p = node->start; p != node->end; ++p) {
                if ((*p)->has_intersection(ray)) {
                    return true;
                }
            }
            return false;
        }
        else
        {
            return has_intersection(ray, node->l) || has_intersection(ray, node->r);
        }
    }
    else
    {
        return false;
    }
}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
    // TODO (Part 2.3):
    // Fill in the intersect function.

    if (node->bb.intersect(ray, ray.min_t, ray.max_t)) {

        // Node is a leaf
        if (node->isLeaf()) {

            bool hit = false;
            for (auto p = node->start; p != node->end; p++) {
                bool p_hit = (*p)->intersect(ray, i);
                hit = hit || p_hit;
            }
            //std::ofstream file("intersect_log.txt");
            //if (file.is_open()) { 
            //    file << i->t << std::endl; 
            //    file.close(); // ¹Ø±ÕÎÄ¼þ
            //}
            return hit;
        }
        else {
            bool hit_left = intersect(ray, i, node->l);
            bool hit_right = intersect(ray, i, node->r);
            return hit_left || hit_right;
        }

    }
    else {
        return false;
    }

}

} // namespace SceneObjects
} // namespace CGL
