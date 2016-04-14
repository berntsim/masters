#include "containers.h"
#include "routines.h"
#include <iostream>

bool AABB::contains(Point &a) const {                                           //This routine checks if point a is located within a certain axis aligned bounded
    if ((a.x < bottom_right.x) && (a.x > top_left.x)){      //box. Here the coordinates are compared to see if the x-component of a is within
        if ((a.y < top_left.y) && (a.y > bottom_right.y)){  //the regioned spanned by the AABB. Silimarly, one checks the y-direction.
            return true;                                                        //If the point is within both regions, it is within the box and we return true.
        }
    }
    return false;                                                               // If point a is not within the box, we return false.
}

bool AABB::intersects(AABB &other) const {
    if ((top_left.x <= other.bottom_right.x) &&
        (bottom_right.x >= other.top_left.x) &&
        (top_left.y >= other.bottom_right.y) &&
        (bottom_right.y <= other.top_left.y)){
        return true;
    }
    else {
        return false;
    }
}

//bool AABB::intersects(AABB &other) const {
//    if (((this->bottom_right.x < this->top_left.x) ||
//        (this->bottom_right.y > this->top_left.y)) ||
//        ((other.bottom_right.x < other.top_left.x) ||
//         (other.bottom_right.y > other.top_left.y))){
//        Point tl = this->top_left;
//        Point br = this->bottom_right;
//        Point tr = Point(this->bottom_right.x, this->top_left.y);
//        Point bl = Point(this->top_left.x, this->bottom_right.y);
//        if ((checkPointInAABB(tl, other))){
//            return true;
//        }
//        else if (checkPointInAABB(tr, other)){
//            return true;
//        }
//        else if (checkPointInAABB(bl, other)){
//            return true;
//        }
//        else if (checkPointInAABB(br, other)){
//            return true;
//        }
//        else if (checkPointInAABB(other.top_left, *this)){
//            return true;
//        }
//        else if (checkPointInAABB(other.bottom_right, *this)){
//            return true;
//        }
//        else {
//            return false;
//        }
//    }
//    else {
//        if ((top_left.x < other.bottom_right.x) &&
//            (bottom_right.x > other.top_left.x) &&
//            (top_left.y > other.bottom_right.y) &&
//            (bottom_right.y < other.top_left.y)){
//            return true;
//        }
//        else {
//            return false;
//        }
//    }
//}

Quadtree::Quadtree(AABB _boundary) :
    nw{nullptr},
    ne{nullptr},
    sw{nullptr},
    se{nullptr},
    boundary(_boundary){
    objects = std::vector<Cluster>();
}

Quadtree::~Quadtree(){
    delete nw;
    delete sw;
    delete ne;
    delete se;
}

//void Quadtree::subdivide(){
//    Point tl = boundary.top_left;
//    Point br = Point((boundary.bottom_right.x - boundary.top_left.x)/2.0,
//                     (boundary.top_left.y - boundary.top_left.y)/2.0);
//    nw = new Quadtree(AABB(tl, br));

//    tl = br;
//    br = boundary.bottom_right;
//    se = new Quadtree(AABB(tl, br));

//    tl.x = (boundary.bottom_right.x - boundary.top_left.x)/2.0;
//    tl.y = boundary.top_left.y;
//    br.x = boundary.bottom_right.x;
//    br.y = (boundary.top_left.y - boundary.bottom_right.y)/2.0;
//    ne = new Quadtree(AABB(tl, br));

//    tl.x = boundary.top_left.x;
//    tl.y = (boundary.top_left.y - boundary.bottom_right.y)/2.0;
//    br.x = (boundary.bottom_right.x - boundary.top_left.x);
//    br.y = boundary.bottom_right.y;
//    sw = new Quadtree(AABB(tl, br));
//}

void Quadtree::subdivide(){
    Point tl = boundary.top_left;
    Point br = Point(boundary.top_left.x +
                     (boundary.bottom_right.x - boundary.top_left.x)/2.0,
                     boundary.bottom_right.y +
                     (boundary.top_left.y - boundary.bottom_right.y)/2.0);
    nw = new Quadtree(AABB(tl, br));

    tl.x = boundary.top_left.x +
            (boundary.bottom_right.x - boundary.top_left.x)/2.0;
    tl.y = boundary.bottom_right.y +
            (boundary.top_left.y - boundary.bottom_right.y)/2.0;
    br = boundary.bottom_right;
    se = new Quadtree(AABB(tl, br));

    tl.x = boundary.top_left.x + (boundary.bottom_right.x - boundary.top_left.x)/2.0;
    tl.y = boundary.top_left.y;
    br.x = boundary.bottom_right.x;
    br.y = boundary.bottom_right.y + (boundary.top_left.y - boundary.bottom_right.y)/2.0;
    ne = new Quadtree(AABB(tl, br));

    tl.x = boundary.top_left.x;
    tl.y = boundary.bottom_right.y + (boundary.top_left.y - boundary.bottom_right.y)/2.0;
    br.x = boundary.top_left.x + (boundary.bottom_right.x - boundary.top_left.x)/2.0;
    br.y = boundary.bottom_right.y;
    sw = new Quadtree(AABB(tl, br));
}

bool Quadtree::insert(Cluster d){
    for (auto&& area:d.areas){
        if (!boundary.intersects(area)){
            return false;
        }
    }
    if (objects.size() < capacity){
        objects.push_back(d);
        return true;
    }
    if (nw == nullptr){
        subdivide();
    }
    if (nw->insert(d)){
        return true;
    }
    if (ne->insert(d)){
        return true;
    }
    if (sw->insert(d)){
        return true;
    }
    if (se->insert(d)){
        return true;
    }
    return false;
}

std::vector<Cluster> Quadtree::queryRange(AABB range){                          //This routine takes in a region defined by AABB and returns any points within
    std::vector<Cluster> pInRange = std::vector<Cluster>();                     //that region. These points are stored in a vector.

    if(!boundary.intersects(range)){                                            //If the region we search is outside the domain, the function return the empty
        return pInRange;                                                        //vector of points.
    }

    for(auto&& object: objects){
        for (auto&& area:object.areas){
            if(range.intersects(area)){
                pInRange.push_back(object);
            }
        }
    }

    if(nw == nullptr){                                                          //Now we check if the node is a leaf node (could chekc any child node really).
        return pInRange;                                                        //If it is a leaf node, then we return the list of objects.
   }

    std::vector<Cluster> temp = nw->queryRange(range);                          //If the node is not a leaf node, then one must add the objects from the children
    pInRange.insert(pInRange.end(), temp.begin(), temp.end());                  //This is done by calling the queryRange recursively, until all nodes have been
                                                                                //checked.
    temp = ne->queryRange(range);
    pInRange.insert(pInRange.end(), temp.begin(), temp.end());

    temp = sw->queryRange(range);
    pInRange.insert(pInRange.end(), temp.begin(), temp.end());

    temp = se->queryRange(range);
    pInRange.insert(pInRange.end(), temp.begin(), temp.end());

    return pInRange;                                                            //Finally, the complete list of objects needed for checking is returned.
}

void Quadtree::clearQadtree(){
    if (nw == nullptr){
        objects.clear();
        delete nw;
    }
    else {
        objects.clear();
        nw->clearQadtree();
    }
    if (ne == nullptr){
        objects.clear();
        delete ne;
    }
    else {
        objects.clear();
        ne->clearQadtree();
    }
    if (sw == nullptr){
        objects.clear();
        delete sw;
    }
    else {
        objects.clear();
        sw->clearQadtree();
    }
    if (se == nullptr){
        objects.clear();
        delete se;
    }
    else {
        objects.clear();
        se->clearQadtree();
    }
}
