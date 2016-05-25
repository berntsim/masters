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

bool AABB_vel::contains(Point &a) const {
    if ((a.x < bottom_right.x) && (a.x > top_left.x)){
        if ((a.y < top_left.y) && (a.y > bottom_right.y)){
            return true;
        }
    }
    return false;
}

Quadtree_vel::Quadtree_vel(AABB_vel _boundary, eddy _E, int _level,
                           double _lifetime, double _turnovertime,
                           std::vector<double> _p) :
    nw{nullptr},
    ne{nullptr},
    sw{nullptr},
    se{nullptr},
    boundary(_boundary){
    E = _E;
    level = _level;
    lifetime = _lifetime;
    turnovertime = _turnovertime;
    p = _p;
}

Quadtree_vel::~Quadtree_vel(){
    delete nw;
    delete sw;
    delete ne;
    delete se;
}

void Quadtree_vel::subdivide_vel(std::mt19937 g, std::mt19937::result_type seed,
                                 double PI, int max_lvl, double L_typical){
    std::vector<double> list_copy = this->p;                                        //list is a list of 4 elements, 2 of p_1 and 2 of p_2
    int i_max = std::numeric_limits<int>::max();
    auto rand_dir = std::bind(std::uniform_real_distribution<double>(0,2*PI),
                   std::mt19937(seed));
    auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
                       std::mt19937(seed));
    double sf = 1.0;                                                            //Share factor. The article states "shared equally", suggesting that there should
    if (this->level == max_lvl){                                                //be a share factor of 0.5. This would however lead to an increase in lifetimes
        return;                                                                 //for smaller eddies, which is unphysical.
    }
    std::shuffle(this->p.begin(), this->p.end(), g);
    Point tl = boundary.top_left;
    Point br = Point(boundary.top_left.x +
                     (boundary.bottom_right.x - boundary.top_left.x)/2.0,
                     boundary.bottom_right.y +
                     (boundary.top_left.y - boundary.bottom_right.y)/2.0);
    eddy tmp_eddy = eddy(this->p[this->p.size()-1]*sf*this->E.E, rand_dir(),
                         this->E.L/2.0);
    nw = new Quadtree_vel(AABB_vel(tl, br), tmp_eddy, this->level+1, 0.0,
                          findLifetime(tmp_eddy, L_typical), list_copy);
    nw->subdivide_vel(g, rand_seed(), PI, max_lvl, L_typical);
    this->p.pop_back();

    tl.x = boundary.top_left.x +
            (boundary.bottom_right.x - boundary.top_left.x)/2.0;
    tl.y = boundary.bottom_right.y +
            (boundary.top_left.y - boundary.bottom_right.y)/2.0;
    br = boundary.bottom_right;
    tmp_eddy = eddy(this->p[this->p.size()-1]*sf*this->E.E, rand_dir(),
                    this->E.L/2.0);
    se = new Quadtree_vel(AABB_vel(tl, br), tmp_eddy, this->level+1, 0.0,
                          findLifetime(tmp_eddy, L_typical), list_copy);
    se->subdivide_vel(g, rand_seed(), PI, max_lvl, L_typical);
    this->p.pop_back();

    tl.x = boundary.top_left.x +
           (boundary.bottom_right.x - boundary.top_left.x)/2.0;
    tl.y = boundary.top_left.y;
    br.x = boundary.bottom_right.x;
    br.y = boundary.bottom_right.y +
           (boundary.top_left.y - boundary.bottom_right.y)/2.0;
    tmp_eddy = eddy(this->p[this->p.size()-1]*sf*this->E.E, rand_dir(),
                    this->E.L/2.0);
    ne = new Quadtree_vel(AABB_vel(tl, br), tmp_eddy, this->level+1, 0.0,
                          findLifetime(tmp_eddy, L_typical), list_copy);
    ne->subdivide_vel(g, rand_seed(), PI, max_lvl, L_typical);
    this->p.pop_back();

    tl.x = boundary.top_left.x;
    tl.y = boundary.bottom_right.y +
           (boundary.top_left.y - boundary.bottom_right.y)/2.0;
    br.x = boundary.top_left.x +
           (boundary.bottom_right.x - boundary.top_left.x)/2.0;
    br.y = boundary.bottom_right.y;
    tmp_eddy = eddy(this->p[this->p.size()-1]*sf*this->E.E, rand_dir(),
                    this->E.L/2.0);
    sw = new Quadtree_vel(AABB_vel(tl, br), tmp_eddy, this->level+1, 0.0,
                          findLifetime(tmp_eddy, L_typical), list_copy);
    sw->subdivide_vel(g, rand_seed(), PI, max_lvl, L_typical);
    this->p.pop_back();
}

Velocity Quadtree_vel::queryRange_vel(Point p){                                 //This routine takes in a point, and returns the velocity in that point.
    Velocity V_ret;                                                             //The value is initially zero.
    double vx = 0;
    double vy = 0;

    if(!boundary.contains(p)){                                                  //If the point is outside the domain spanned by the quadtree, then it
        return V_ret;                                                           //returns the velocity vector.
    }
    else {
        vx += std::sqrt(this->E.E)*std::cos(this->E.theta);                     //The contribution from this eddy is added to the total velocity
        vy += std::sqrt(this->E.E)*std::sin(this->E.theta);
        V_ret.v = std::sqrt(vx*vx + vy*vy);
        V_ret.theta = std::tan(vy/vx);
    }

    if(nw == nullptr){                                                          //Now we check if the node is a leaf node (could check any child node really).
        return V_ret;                                                           //If it is a leaf node, then we return the velocity.
   }


    Velocity V_tmp = nw->queryRange_vel(p);                                     //If the node is not a leaf node, then one must add the eddies from the children
    vx += V_tmp.v*std::cos(V_tmp.theta);
    vy += V_tmp.v*std::sin(V_tmp.theta);                                        //This is done by calling the queryRange recursively, until all nodes have been
                                                                                //checked.
    V_tmp = ne->queryRange_vel(p);
    vx += V_tmp.v*std::cos(V_tmp.theta);
    vy += V_tmp.v*std::sin(V_tmp.theta);

    V_tmp = sw->queryRange_vel(p);
    vx += V_tmp.v*std::cos(V_tmp.theta);
    vy += V_tmp.v*std::sin(V_tmp.theta);

    V_tmp = se->queryRange_vel(p);
    vx += V_tmp.v*std::cos(V_tmp.theta);
    vy += V_tmp.v*std::sin(V_tmp.theta);

    V_ret.v = std::sqrt(vx*vx + vy*vy);
    if (std::abs(vx) < std::abs(vy)){
        V_ret.theta = std::asin(vy/V_ret.v);
    }
    else {
        V_ret.theta = std::acos(vx/V_ret.v);
    }
    return V_ret;
}

void Quadtree_vel::clearQuadtree_vel(){
    if (nw == nullptr){
        delete nw;
    }
    else {
        nw->clearQuadtree_vel();
    }
    if (ne == nullptr){
        delete ne;
    }
    else {
        ne->clearQuadtree_vel();
    }
    if (sw == nullptr){
        delete sw;
    }
    else {
        sw->clearQuadtree_vel();
    }
    if (se == nullptr){
        delete se;
    }
    else {
        se->clearQuadtree_vel();
    }
}


void Quadtree_vel::updateEddyTree(double time, int max_lvl, std::mt19937 g,
                                  double PI, std::mt19937::result_type seed,
                                  std::vector<double> p_list, double L_typical){
    if (this->level < max_lvl){
        std::vector<double> list_copy = p_list;
        int i_max = std::numeric_limits<int>::max();
        auto rand_dir = std::bind(std::uniform_real_distribution<double>(0,2*PI),
                       std::mt19937(seed));
        auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
                           std::mt19937(seed));
//        double sf = 1.0;
        if ((time - nw->lifetime) > nw->turnovertime){
            if (int(this->p.size()) > 0){
//                nw->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                nw->E.theta = rand_dir();
                nw->lifetime = time;
                nw->turnovertime = findLifetime(nw->E, L_typical);
                this->p.pop_back();
            }
            else {
                std::shuffle(p_list.begin(), p_list.end(), g);
                this->p = p_list;
//                nw->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                nw->E.theta = rand_dir();
                nw->lifetime = time;
                nw->turnovertime = findLifetime(nw->E, L_typical);
                this->p.pop_back();
            }
        }
        if ((time - ne->lifetime) > ne->turnovertime){
            if (int(this->p.size()) > 0){
//                ne->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                ne->E.theta = rand_dir();
                ne->lifetime = time;
                ne->turnovertime = findLifetime(ne->E, L_typical);
                this->p.pop_back();
            }
            else {
                std::shuffle(p_list.begin(), p_list.end(), g);
                this->p = p_list;
//                ne->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                ne->E.theta = rand_dir();
                ne->lifetime = time;
                ne->turnovertime = findLifetime(ne->E, L_typical);
                this->p.pop_back();
            }
        }
        if ((time - sw->lifetime) > sw->turnovertime){
            if (int(this->p.size()) > 0){
//                sw->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                sw->E.theta = rand_dir();
                sw->lifetime = time;
                sw->turnovertime = findLifetime(sw->E, L_typical);
                this->p.pop_back();
            }
            else {
                std::shuffle(p_list.begin(), p_list.end(), g);
                this->p = p_list;
//                sw->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                sw->E.theta = rand_dir();
                sw->lifetime = time;
                sw->turnovertime = findLifetime(sw->E, L_typical);
                this->p.pop_back();
            }
        }
        if ((time - se->lifetime) > se->turnovertime){
            if (int(this->p.size()) > 0){
//                se->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                se->E.theta = rand_dir();
                se->lifetime = time;
                se->turnovertime = findLifetime(se->E, L_typical);
                this->p.pop_back();
            }
            else {
                std::shuffle(p_list.begin(), p_list.end(), g);
                this->p = p_list;
//                se->E.E = this->p[this->p.size()-1]*sf*this->E.E;
                se->E.theta = rand_dir();
                se->lifetime = time;
                se->turnovertime = findLifetime(se->E, L_typical);
                this->p.pop_back();
            }
        }
        nw->updateEddyTree(time, max_lvl, g, PI, rand_seed(), list_copy, L_typical);
        ne->updateEddyTree(time, max_lvl, g, PI, rand_seed(), list_copy, L_typical);
        sw->updateEddyTree(time, max_lvl, g, PI, rand_seed(), list_copy, L_typical);
        se->updateEddyTree(time, max_lvl, g, PI, rand_seed(), list_copy, L_typical);
    }
}

