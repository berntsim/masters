#ifndef VEL_FIELD_H
#define VEL_FIELD_H

#include <vector>
#include <random>
#include <cmath>
#include "containers.h"

struct eddy{
    double E, theta, L;
    eddy(double _E = 0, double _theta = 0, double _L = 0) :
        E(_E), theta(_theta), L(_L){}
};

struct AABB_vel{
    Point top_left;
    Point bottom_right;
    AABB_vel(Point _top_left, Point _bottom_right) : top_left(_top_left),
             bottom_right(_bottom_right){}
    bool contains(Point &a) const;
};

class Quadtree_vel{
    public:
        Quadtree_vel* nw;                                                           //Pointer to the Quadtree node for the nw quarter partitioning of the region.
        Quadtree_vel* ne;                                                           //Similar, only the ne quarter
        Quadtree_vel* sw;                                                           //Similar, only the sw quarter
        Quadtree_vel* se;                                                           //Similar, only the se quarter
        AABB_vel boundary;                                                          //Boundary of the space partitioned in this node.
//        velocity E;
        eddy E;
        int level;
        double lifetime;
        double turnovertime;
        std::vector<float> p;

//------Methods-----------------------------------------------------------------//This is what would usually be public
        Quadtree_vel();
        Quadtree_vel(AABB_vel _boundary, eddy E, int _level,
                     double _lifetime, double _turnovertime,
                     std::vector<float> _p);

        ~Quadtree_vel();

        void subdivide_vel(std::mt19937 g, std::mt19937::result_type seed,
                           double PI, int max_lvl);
        velocity queryRange_vel(Point p);
        void clearQuadtree_vel();
        void updateEddyTree(double time, int max_lvl, std::mt19937 g,
                            double PI, std::mt19937::result_type seed,
                            std::vector<float> p_list);
};

#endif // VEL_FIELD_H
