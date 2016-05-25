#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>

struct Point{
    double x, y;
    Point(double _x = 0, double _y = 0) : x(_x), y(_y){}
};

struct AABB{
    Point top_left;
    Point bottom_right;
    AABB(Point _top_left, Point _bottom_right) : top_left(_top_left),
        bottom_right(_bottom_right){}
    bool contains(Point &a) const;
    bool intersects(AABB &other) const;
};

struct Particle{
    Point pos;
    double r_p;
    Particle(Point _pos, double _r_p): pos(_pos), r_p(_r_p){}

};

struct Velocity{
    double v, theta;
    Velocity(double _v = 0, double _theta = 0) : v(_v), theta(_theta){}
};

struct Cluster{
    bool joined;
    std::vector<AABB> areas;
    std::vector<Particle> particles;
    int index;
    double radius;
    Velocity vel;
    Point CM;
    double mass;
    double lifetime;
    Cluster(bool _joined, std::vector<AABB> _areas, std::vector<Particle> _particles,
            int _index, double _radius, Velocity _vel, Point _CM, double _mass,
            double _lifetime):
        joined(_joined), areas(_areas), particles(_particles), index(_index),
        radius(_radius), vel(_vel), CM(_CM), mass(_mass), lifetime(_lifetime){}
};

class Quadtree{
    private:
        Quadtree* nw;                                                           //Pointer to the Quadtree node for the nw quarter partitioning of the region.
        Quadtree* ne;                                                           //Similar, only the ne quarter
        Quadtree* sw;                                                           //Similar, only the sw quarter
        Quadtree* se;                                                           //Similar, only the se quarter
        AABB boundary;                                                          //Boundary of the space partitioned in this node.


    public:
        std::vector<Cluster> objects;                                           //These objects should be particles in my case. Note that it would probably make
                                                                                //sense to not have this as a vector, since I have capacity = 1.
        static constexpr int capacity = 225;                                    //Setting that every node can at contain at MAXIMUM 1 particle.
        Quadtree();                                                             //Default constructor. This will set all children to 0, default AABB and objects.
        Quadtree(AABB _boundary);                                               //With this constructor, the boundary is determined by the argument.

        ~Quadtree();                                                            //Destructor deleting all the dynamically allocated pointers (nw,ne,sw,se).

        bool insert(Cluster d);                                                 //Will insert a cluster object into the appropriate place in the quadtree.
        void subdivide();                                                       //Will split the quadtree to four new regions, making the quadtree recursive.
        std::vector<Cluster> queryRange(AABB range);                            //A function searching the given range for all cluster ovjects who overlap with
        void clearQadtree();                                                    //this area. These are then put into a list and returned.

};

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
        std::vector<double> p;

//------Methods-----------------------------------------------------------------//This is what would usually be public
        Quadtree_vel();
        Quadtree_vel(AABB_vel _boundary, eddy E, int _level,
                     double _lifetime, double _turnovertime,
                     std::vector<double> _p);

        ~Quadtree_vel();

        void subdivide_vel(std::mt19937 g, std::mt19937::result_type seed,
                           double PI, int max_lvl, double L_typical);
        Velocity queryRange_vel(Point p);
        void clearQuadtree_vel();
        void updateEddyTree(double time, int max_lvl, std::mt19937 g,
                            double PI, std::mt19937::result_type seed,
                            std::vector<double> p_list, double L_typical);
};




#endif // CONTAINERS_H

