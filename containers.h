#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <vector>
#include <iostream>

struct Point{
    float x, y;
    Point(float _x = 0, float _y = 0) : x(_x), y(_y){}
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

struct Cluster{
    bool joined;
    AABB area;
    int index;
    std::vector<Particle> particles;
    Cluster(bool _joined, AABB _area, std::vector<Particle> _particles,
            int _index):
        joined(_joined), area(_area), particles(_particles), index(_index){}
};

class Quadtree{
    private:
        Quadtree* nw;                                                           //Pointer to the Quadtree node for the nw quarter partitioning of the region.
        Quadtree* ne;                                                           //Similar, only the ne quarter
        Quadtree* sw;                                                           //Similar, only the sw quarter
        Quadtree* se;                                                           //Similar, only the se quarter
        AABB boundary;                                                          //Boundary of the space partitioned in this node.
        std::vector<Cluster> objects;                                           //These objects should be particles in my case. Note that it would probably make
                                                                                //sense to not have this as a vector, since I have capacity = 1.
        static constexpr int capacity = 225;                                      //Setting that every node can at contain at MAXIMUM 1 particle.
    public:
        Quadtree();                                                             //Default constructor. This will set all children to 0, default AABB and objects.
        Quadtree(AABB _boundary);                                               //With this constructor, the boundary is determined by the argument.

        ~Quadtree();                                                            //Destructor deleting all the dynamically allocated pointers (nw,ne,sw,se).

        bool insert(Cluster d);                                                 //Will insert a cluster object into the appropriate place in the quadtree.
        void subdivide();                                                       //Will split the quadtree to four new regions, making the quadtree recursive.
        std::vector<Cluster> queryRange(AABB range);                            //A function searching the given range for all cluster ovjects who overlap with
        void clearQadtree();                                                           //this area. These are then put into a list and returned.

};




#endif // CONTAINERS_H