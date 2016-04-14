#include "routines.h"
#include "containers.h"

double findDistance(Point a, Point b, double x_size, double y_size){
    double dx = a.x - b.x;
    if (dx > x_size * 0.5){
        dx = dx - x_size;
    }
    else if (dx <= -x_size * 0.5){
        dx = dx + x_size;
    }

    double dy = a.y - b.y;

    if (dy > y_size * 0.5){
        dy = dy - y_size;
    }
    else if (dy <= -y_size * 0.5){
        dy = dy + y_size;
    }
    return std::sqrt(pow(dx,2) + pow(dy,2));
}


void distributeParticlesTest(int len, int height, int nbr_particles,
                            std::vector<Cluster> &clusters,
                            std::mt19937::result_type seed_x,
                            std::mt19937::result_type seed_y, double r_p,
                            bool varying_size, std::map<int, Cluster*> &clust_test){
    auto coord_x = std::bind(std::uniform_real_distribution<double>(0,len),     //We define a function to generate a random point in the x-dimension on the domain
                               std::mt19937(seed_x));
    auto coord_y = std::bind(std::uniform_real_distribution<double>(0,height),  //similarly for the y-dimension.
                               std::mt19937(seed_y));
    auto rand_size = std::bind(std::uniform_real_distribution<float>(1.0,5.0),
                               std::mt19937(seed_x));
    double dist;
    int part_placed = 0;
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;


    while (int(clusters.size()) < nbr_particles){
        tmp.x = coord_x();
        tmp.y = coord_y();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clusters.size()); ++i){
            dist = findDistance(clusters[i].particles[0].pos, tmp, len, height);
            if (clusters[i].particles[0].r_p + r_p > dist){
                occupied = true;
            }
        }
        if (!occupied){
            tmp_part.push_back(Particle(tmp, r_p));
            aabb_tmp.bottom_right.x = tmp.x+r_p;
            aabb_tmp.bottom_right.y = tmp.y - r_p;
            aabb_tmp.top_left.x = tmp.x - r_p;
            aabb_tmp.top_left.y = tmp.y + r_p;
            areas_tmp.push_back(aabb_tmp);
//            Cluster tmp_clust = Cluster(aabb_tmp, tmp_part, part_placed);
            clusters.push_back(Cluster(false, areas_tmp, tmp_part, part_placed));
            clust_test.insert(std::make_pair(part_placed, new Cluster(false, areas_tmp, tmp_part, part_placed)));
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;                                                       //We reset and prepare to place the next particle.
    }
}




void addToDraw(std::vector<sf::CircleShape> &to_draw,
               std::vector<Cluster> clusters,
               sf::Color fill_color){
    sf::CircleShape shape(1.0);
    for (auto&& clust: clusters){
        for (auto&& particle: clust.particles){
            shape.setRadius(particle.r_p);
            shape.setPosition(particle.pos.x - particle.r_p,
                              particle.pos.y - particle.r_p);
            shape.setFillColor(fill_color);
            to_draw.push_back(shape);
        }
    }
}



sf::VertexArray printSearchRange(AABB range){
    sf::VertexArray lines(sf::LinesStrip, 5);
    lines[0].position = sf::Vector2f(range.top_left.x, range.top_left.y);
    lines[1].position = sf::Vector2f(range.bottom_right.x, range.top_left.y);
    lines[2].position = sf::Vector2f(range.bottom_right.x,range.bottom_right.y);
    lines[3].position = sf::Vector2f(range.top_left.x, range.bottom_right.y);
    lines[4].position = sf::Vector2f(range.top_left.x, range.top_left.y);
    return lines;
}

void putInQuadtree(std::vector<Cluster> clusters, Quadtree &qtree){
    for (auto&& cluster: clusters){
        qtree.insert(cluster);
    }
}

void testPutInQuadtree(std::map<int, Cluster*> clusters, Quadtree &qtree){
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end(); iterator++){
        if (iterator->second != nullptr){
            qtree.insert(*(iterator->second));
        }
    }
}

void takeStep(std::mt19937::result_type seed, double len,
              std::vector<Cluster> &clusters,int x_size, int y_size){
    double PI = 3.1415926535897932384626433832;
    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*PI),
                   std::mt19937(seed));
    double dir;
    for (auto&& cluster: clusters){
        dir = rand_dir();
        for (auto&& particle: cluster.particles){
            particle.pos.x += len*std::cos(dir);
            particle.pos.y += len*std::sin(dir);
            if (particle.pos.x < 0){
                particle.pos.x += x_size;
            }
            else if (particle.pos.x >= x_size){
                particle.pos.x -= x_size;
            }
            if (particle.pos.y < 0){
                particle.pos.y += y_size;
            }
            else if (particle.pos.y > y_size){
                particle.pos.y -= y_size;
            }
        }
        cluster.areas[0].top_left.x += len*std::cos(dir);
        cluster.areas[0].bottom_right.x += len*std::cos(dir);
        cluster.areas[0].top_left.y += len*std::sin(dir);
        cluster.areas[0].bottom_right.y += len*std::sin(dir);
    }
}


void takeSingleStep(double step_dir, double len, Cluster* cluster, int x_size,
                    int y_size){
    if (cluster != nullptr){
        for (auto&& particle: cluster->particles){
            particle.pos.x += len*std::cos(step_dir);
            particle.pos.y += len*std::sin(step_dir);
            if (particle.pos.x < 0){
                particle.pos.x += x_size;
            }
            else if (particle.pos.x >= x_size){
                particle.pos.x -= x_size;
            }
            if (particle.pos.y < 0){
                particle.pos.y += y_size;
            }
            else if (particle.pos.y >= y_size){
                particle.pos.y -= y_size;
            }
        }
        for (auto&& area:cluster->areas){
            area.top_left.x += len*std::cos(step_dir);
            area.bottom_right.x += len*std::cos(step_dir);
            area.top_left.y += len*std::sin(step_dir);
            area.bottom_right.y += len*std::sin(step_dir);
            if (area.top_left.x < 0){
                area.top_left.x += x_size;
            }
            else if (area.top_left.x > x_size){
                area.top_left.x -= x_size;
            }
            if (area.bottom_right.x < 0){
                area.bottom_right.x += x_size;
            }
            else if (area.bottom_right.x > x_size){
                area.bottom_right.x -= x_size;
            }
            if (area.top_left.y > y_size){
                area.top_left.y -= y_size;
            }
            else if (area.top_left.y < 0){
                area.top_left.y += y_size;
            }
            if (area.bottom_right.y > y_size){
                area.bottom_right.y -= y_size;
            }
            else if (area.bottom_right.y < 0){
                area.bottom_right.y += y_size;
            }
        }
    }
}

void findSearchRange(std::vector<AABB> &search_ranges, Cluster* cluster,
                     double step_len, double x_size, double y_size){
    if (cluster != nullptr){
        AABB range = AABB(Point(0,0), Point(0,0));
        double tlx, tly, brx, bry;
        tlx = cluster->areas[0].top_left.x - step_len;
        tly = cluster->areas[0].top_left.y + step_len;
        brx = cluster->areas[0].bottom_right.x + step_len;
        bry = cluster->areas[0].bottom_right.y - step_len;

        if (tlx < 0){
            tlx += x_size;
            if (tly > y_size){
                tly -= y_size;
                range.top_left.x = tlx;
                range.top_left.y = tly;
                range.bottom_right.x = x_size;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (1 in notes)
                range.top_left.y = y_size;
                range.bottom_right.x = x_size;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (2 in notes)
                range.top_left.x = 0;
                range.top_left.y = tly;
                range.bottom_right.x = brx;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (6 in notes)
                range.top_left.y = y_size;
                range.bottom_right.x = brx;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (7 in notes)
            }
            else if (bry < 0){
                bry += y_size;
                range.top_left.x = tlx;
                range.top_left.y = y_size;
                range.bottom_right.x = x_size;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (3 in notes)
                range.top_left.x = tlx;
                range.top_left.y = tly;
                range.bottom_right.x = x_size;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (4 in notes)
                range.top_left.x = 0;
                range.top_left.y = y_size;
                range.bottom_right.x = brx;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (8 in notes)
                range.top_left.y = tly;
                range.bottom_right.x = brx;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (9 in notes)
            }
            else {
                range.top_left.x = tlx;
                range.top_left.y = tly;
                range.bottom_right.x = x_size;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (5 in notes)
                range.top_left.x = 0;
                range.bottom_right.x = brx;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (10 in notes)
            }
        }
        else if (brx >= x_size){
            brx -= x_size;
            if (tly > y_size){
                tly -= y_size;
                range.top_left.x = 0;
                range.top_left.y = tly;
                range.bottom_right.x = brx;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (11 in notes)
                range.top_left.y = y_size;
                range.bottom_right.x = brx;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (12 in notes)
                range.top_left.x = tlx;
                range.top_left.y = tly;
                range.bottom_right.x = x_size;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (16 in notes)
                range.top_left.y = y_size;
                range.bottom_right.x = x_size;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (17 in notes)
            }
            else if (bry < 0){
                bry += y_size;
                range.top_left.x = 0;
                range.top_left.y = y_size;
                range.bottom_right.x = brx;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (13 in notes)
                range.top_left.x = 0;
                range.top_left.y = tly;
                range.bottom_right.x = brx;
                range.bottom_right.y = 0;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (14 in notes)
                range.top_left.x = tlx;
                range.top_left.y = y_size;
                range.bottom_right.x = x_size;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (18 in notes)
                range.top_left.y = tly;
                range.bottom_right.x = x_size;
                range.bottom_right.y = y_size;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (19 in notes)
            }
            else {
                range.top_left.x = 0;
                range.top_left.y = tly;
                range.bottom_right.x = brx;
                range.bottom_right.y = bry;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (15 in notes)
                range.top_left.x = tlx;
                range.bottom_right.x = x_size;
                search_ranges.push_back(range);                                 //Insert search range into search_ranges (20 in notes)
            }

        }
        else if (bry < 0){
            bry += y_size;
            range.top_left.x = tlx;
            range.top_left.y = y_size;
            range.bottom_right.x = brx;
            range.bottom_right.y = bry;
            search_ranges.push_back(range);                                     //Insert search range into search_ranges (21 in notes)
            range.top_left.y = tly;
            range.bottom_right.y = 0;
            search_ranges.push_back(range);                                     //Insert search range into search_ranges (22 in notes)

        }
        else if (tly >= y_size){
           tly -= y_size;
           range.top_left = tly;
           range.top_left.x = tlx;
           range.bottom_right.x = brx;
           range.bottom_right.y = 0;
           search_ranges.push_back(range);                                      //Insert search range into search_ranges (23 in notes)
           range.top_left.y = y_size;
           range.bottom_right.y = bry;
           search_ranges.push_back(range);                                      //Insert search range into search_ranges (24 in notes)
        }
        else {
            range.top_left.x = tlx;
            range.top_left.y = tly;
            range.bottom_right.x = brx;
            range.bottom_right.y = bry;
            search_ranges.push_back(range);
        }
    }
}

std::vector<Cluster> BPcolCheck(Cluster* cluster,
                                std::vector<AABB> search_ranges,
                                Quadtree &tree){
    std::vector<Cluster> ret_vec;
    if (cluster != nullptr){
        std::vector<Cluster> tmp_res;
        int index = cluster->index;
        for (auto&& search_range:search_ranges){
            tmp_res = tree.queryRange(search_range);
            if (tmp_res.size() > 1){
                for (auto&& res: tmp_res){
                    if (res.index != index){
                        ret_vec.push_back(res);
                    }
                }
            }
        }
    }
    return ret_vec;
}

double LHit(double step_L, double step_dir, Particle one, Particle two,
            int x_size, int y_size){
    //rewrite the check to take findDistance into account!
    double dx = one.pos.x - two.pos.x;
    if (dx > x_size * 0.5){
        dx = dx - x_size;
    }
    else if (dx <= -x_size * 0.5){
        dx = dx + x_size;
    }

    double dy = one.pos.y - two.pos.y;

    if (dy > y_size * 0.5){
        dy = dy - y_size;
    }
    else if (dy <= -y_size * 0.5){
        dy = dy + y_size;
    }
//    std::cout << std::endl;
//    std::cout << "pos one = " << one.pos.x << "," << one.pos.y << std::endl;
//    std::cout << "pos two = " << two.pos.x << "," << two.pos.y << std::endl;
//    std::cout << "dx = " << dx << ", dy = " << dy << std::endl;
//    std::cout << "step_dir = " << step_dir << std::endl;
    double d_p = one.r_p + two.r_p;
    double a = 1.0;
    double b = 2.0*(std::cos(step_dir)*(dx) +
                    std::sin(step_dir)*(dy));
    double c = (-dx)*(-dx) +
               (-dy)*(-dy) - d_p*d_p;
    double res[2];
    bool sol = true;
//    std::cout << 2*std::sin(-3.1415926535/2.0)*dy << std::endl;
//    std::cout << "a = " << a << std::endl;
//    std::cout << "b = " << b << std::endl;
//    std::cout << "c = " << c << std::endl;
//    std::cout << "d_p = " << d_p << std::endl;
    if ((b*b - 4*a*c) < 0){
//        std::cout << "False solution" << std::endl;
        sol = false;
        return step_L;
    }
    else {
        res[0] = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
        res[1] = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
//        std::cout << "res[0] = " << res[0] << std::endl;
//        std::cout << "res[1] = " << res[1] << std::endl;
        if ((res[0] < res[1]) && ((res[0] > 0) && (res[0] < step_L))){
            return res[0];
        }
        else if ((res[1] < res[0]) && ((res[1] > 0) && (res[1] < step_L))){
            return res[1];
        }
        else {
            return step_L;
        }
    }
}


double NPColCheck(Cluster* cluster, std::vector<Cluster> targets,
                  double step_len, double step_dir, int &col_with,
                  int x_size, int y_size){
    double tmp;
    double ret_val = step_len;
    if (cluster != nullptr){
        for (auto&& target : targets){
            if (target.index != cluster->index){
                for (auto&& particle: cluster->particles){
                    for (auto&& tar_part: target.particles){
                        tmp = LHit(step_len, step_dir, particle, tar_part, x_size,
                                   y_size);
                        if (tmp < ret_val){
                            ret_val = tmp;
                            col_with = target.index;
                        }
                    }
                }
            }
        }
    }
    return ret_val;
}



void joinClusters(Cluster &clust, Cluster &other, std::vector<Cluster> &clusters){
    clust.particles.insert(clust.particles.end(), other.particles.begin(),
                           other.particles.end());
    if (other.areas[0].top_left.x < clust.areas[0].top_left.x){
        clust.areas[0].top_left.y = other.areas[0].top_left.y;
    }
    if (other.areas[0].top_left.y > clust.areas[0].top_left.y){
        clust.areas[0].top_left.y = other.areas[0].top_left.y;
    }
    if (other.areas[0].bottom_right.x > clust.areas[0].bottom_right.x){
        clust.areas[0].bottom_right.x = other.areas[0].bottom_right.x;
    }
    if (other.areas[0].bottom_right.y < clust.areas[0].bottom_right.y){
        clust.areas[0].bottom_right.y = other.areas[0].bottom_right.y;
    }
//    clusters.erase(clusters.begin() + other.index);
}













void testAddToDrawRec(std::vector<sf::RectangleShape> &to_draw,
                      std::map<int, Cluster*> clusters,
                      sf::Color fill_color){
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end(); iterator++){
        sf::RectangleShape shape(sf::Vector2f((iterator->second->areas[0].bottom_right.x - iterator->second->areas[0].top_left.x),
                                              (iterator->second->areas[0].top_left.y - iterator->second->areas[0].bottom_right.y)));
        shape.setPosition(iterator->second->areas[0].top_left.x, iterator->second->areas[0].top_left.y-(iterator->second->areas[0].top_left.y - iterator->second->areas[0].bottom_right.y));
        shape.setFillColor(fill_color);
        to_draw.push_back(shape);
    }
}

void testAddToDraw(std::vector<sf::CircleShape> &to_draw,
               std::map<int, Cluster*> clusters,
               sf::Color fill_color){
    sf::CircleShape shape(1.0);
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end();
        iterator++){
        if (iterator->second != nullptr){
            for (auto&& particle: iterator->second->particles){
                shape.setRadius(particle.r_p);
                shape.setPosition(particle.pos.x - particle.r_p,
                                  particle.pos.y - particle.r_p);
                shape.setFillColor(fill_color);
                to_draw.push_back(shape);
            }
        }
    }
}

void TestJoinClusters(Cluster* clust, Cluster* other,
                      std::map<int, Cluster*> &clusters, int x_size,
                      int y_size){
    if ((clust != nullptr) && (other != nullptr)){
        if (clust->index != other->index){
            updateAreas(clust, other, x_size, y_size);
            clust->particles.insert(clust->particles.end(),
                                    other->particles.begin(),
                                    other->particles.end());
            clusters.at(other->index) = nullptr;
        }
    }
}

void printInformation(Cluster* clust, int col_with){
    for (auto&& particle:clust->particles){
        std::cout << particle.pos.x << ", " << particle.pos.y << "      " << clust->index << "   ";
        std::cout << clust->areas[0].top_left.x << ", " << clust->areas[0].top_left.y << "      ";
        std::cout << clust->areas[0].bottom_right.x << ", " << clust->areas[0].bottom_right.y << "      "
                  << col_with << std::endl;
    }
}


bool crossXBound(Cluster clust, int x_size){
    if (clust.areas.size() == 1){
        if (((clust.areas[0].top_left.x <= x_size) &&
             (clust.areas[0].top_left.x > x_size/2.0)) &&
            ((clust.areas[0].bottom_right.x >= 0) &&
             (clust.areas[0].bottom_right.x < x_size/2.0))){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust.areas.size() == 2){
        if (((clust.areas[0].top_left.x <= x_size) &&
             (clust.areas[0].top_left.x > x_size/2.0)) &&                       // when the separation of the boxes are vertical
            ((clust.areas[1].bottom_right.x >= 0) &&
             (clust.areas[1].bottom_right.x < x_size/2.0))){
            return true;
        }
        else if (((clust.areas[0].top_left.x <= x_size) &&
                  (clust.areas[0].top_left.x > x_size/2.0)) &&                  //Horizontal
                 ((clust.areas[0].bottom_right.x >= 0) &&
                  (clust.areas[0].bottom_right.x < x_size/2.0))){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust.areas.size() == 4){
        if (((clust.areas[0].top_left.x <= x_size) &&
             (clust.areas[0].top_left.x > x_size/2.0)) &&
            ((clust.areas[1].bottom_right.x >= 0) &&
             (clust.areas[1].bottom_right.x < x_size/2.0))){
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

bool crossYBound(Cluster clust, int y_size){
    if (clust.areas.size() == 1){
        if (((clust.areas[0].top_left.y >= 0) &&
             (clust.areas[0].top_left.y < y_size/2.0)) &&
            ((clust.areas[0].bottom_right.y <= y_size) &&
             (clust.areas[0].bottom_right.y > y_size/2.0))){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust.areas.size() == 2){
        if (((clust.areas[0].top_left.y >= 0) &&
             (clust.areas[0].top_left.y < y_size/2.0)) &&                       // when the separation of the boxes are horizontal
            ((clust.areas[1].bottom_right.y <= y_size) &&
             (clust.areas[1].bottom_right.y > y_size/2.0))){
            return true;
        }
        else if (((clust.areas[0].top_left.y >= 0) &&
                  (clust.areas[0].top_left.y < y_size/2.0)) &&                  //vertical
                 ((clust.areas[0].bottom_right.y <= y_size) &&
                  (clust.areas[0].bottom_right.y > y_size/2.0))){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust.areas.size() == 4){
        if (((clust.areas[0].top_left.y >= 0) &&
             (clust.areas[0].top_left.y < y_size/2.0)) &&
            ((clust.areas[2].bottom_right.y <= y_size) &&
             (clust.areas[2].bottom_right.y > y_size/2.0))){
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

void splitAreas(Cluster &clust, int x_size, int y_size){
    if ((crossXBound(clust, x_size)) && !(crossYBound(clust, y_size))){         //The area spanned by the cluster only crosses the x-boundary
        AABB zero(Point(0,0), Point(0,0));
        AABB one(Point(0,0), Point(0,0));
        if (clust.areas.size() == 1){
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right.x = x_size;
            zero.bottom_right.y = clust.areas[0].bottom_right.y;
            one.top_left.x = 0;
            one.top_left.y = clust.areas[0].top_left.y;
            one.bottom_right = clust.areas[0].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
            clust.areas.push_back(one);
        }
        else if (clust.areas.size() == 2){
            if (clust.areas[0].top_left.y == clust.areas[1].top_left.y){
                zero.top_left = clust.areas[0].top_left;
                zero.bottom_right.x = x_size;
                zero.bottom_right.y = clust.areas[0].bottom_right.y;
                one.top_left.x = 0;
                one.top_left.y = clust.areas[0].top_left.y;
                one.bottom_right = clust.areas[1].bottom_right;
                clust.areas.clear();
                clust.areas.push_back(zero);
                clust.areas.push_back(one);
            }
            else {
                zero.top_left = clust.areas[0].top_left;
                zero.bottom_right.x = x_size;
                zero.bottom_right.y = clust.areas[1].bottom_right.y;
                one.top_left.x = 0;
                one.top_left.y = clust.areas[0].top_left.y;
                one.bottom_right = clust.areas[1].bottom_right;
                clust.areas.clear();
                clust.areas.push_back(zero);
                clust.areas.push_back(one);
            }
        }
        else if (clust.areas.size() == 4) {
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right.x = x_size;
            zero.bottom_right.y = clust.areas[3].bottom_right.y;
            one.top_left.x = 0;
            one.top_left.y = clust.areas[1].top_left.y;
            one.bottom_right = clust.areas[2].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
            clust.areas.push_back(one);
        }
    }
    else if ((crossYBound(clust, y_size)) && !(crossXBound(clust, x_size))){    //The area spanned by the cluster only crosses the y-boundary.
        AABB zero(Point(0,0), Point(0,0));
        AABB one(Point(0,0), Point(0,0));
        if (clust.areas.size() == 1){
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right.x = clust.areas[0].bottom_right.x;
            zero.bottom_right.y = 0;
            one.top_left.y = y_size;
            one.top_left.x = clust.areas[0].top_left.x;
            one.bottom_right = clust.areas[0].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
            clust.areas.push_back(one);
        }
        else if (clust.areas.size() == 2){
            if (clust.areas[0].top_left.y == clust.areas[1].top_left.y){
                zero.top_left = clust.areas[0].top_left;
                zero.bottom_right.y = 0;
                zero.bottom_right.x = clust.areas[1].bottom_right.x;
                one.top_left.x = clust.areas[0].top_left.x;
                one.top_left.y = y_size;
                one.bottom_right = clust.areas[1].bottom_right;
                clust.areas.clear();
                clust.areas.push_back(zero);
                clust.areas.push_back(one);
            }
            else {
                zero.top_left = clust.areas[0].top_left;
                zero.bottom_right.x = clust.areas[0].bottom_right.x;
                zero.bottom_right.y = 0;
                one.top_left.y = y_size;
                one.top_left.x = clust.areas[0].top_left.x;
                one.bottom_right = clust.areas[1].bottom_right;
                clust.areas.clear();
                clust.areas.push_back(zero);
                clust.areas.push_back(one);
            }
        }
        else if (clust.areas.size() == 4) {
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right.x = clust.areas[1].bottom_right.x;
            zero.bottom_right.y = 0;
            one.top_left.y = y_size;
            one.top_left.x = clust.areas[3].top_left.x;
            one.bottom_right = clust.areas[2].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
            clust.areas.push_back(one);
        }
    }
    else if ((crossXBound(clust, x_size) && (crossYBound(clust, y_size)))){
        AABB zero(Point(0,0), Point(0,0));
        AABB one(Point(0,0), Point(0,0));
        AABB two(Point(0,0), Point(0,0));
        AABB three(Point(0,0), Point(0,0));
        if (clust.areas.size() == 1){
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right.x = x_size;
            zero.bottom_right.y = 0;
            one.top_left.x = 0;
            one.top_left.y = clust.areas[0].top_left.y;
            one.bottom_right.x = clust.areas[0].bottom_right.x;
            one.bottom_right.y = 0;
            two.top_left.x = 0;
            two.top_left.y = y_size;
            two.bottom_right = clust.areas[0].bottom_right;
            three.top_left.x = clust.areas[0].top_left.x;
            three.top_left.y = y_size;
            three.bottom_right.x = x_size;
            three.bottom_right.y = clust.areas[0].bottom_right.y;
            clust.areas.clear();
            clust.areas.push_back(zero);
            clust.areas.push_back(one);
            clust.areas.push_back(two);
            clust.areas.push_back(three);
        }
        else if (clust.areas.size() == 2){
            if (clust.areas[0].top_left.y == clust.areas[1].top_left.y){        //the areas comming in are split vertically
                zero.top_left = clust.areas[0].top_left;
                zero.bottom_right.x = x_size;
                zero.bottom_right.y = 0;
                one.top_left.x = 0;
                one.top_left.y = clust.areas[1].top_left.y;
                one.bottom_right.x = clust.areas[1].bottom_right.x;
                one.bottom_right.y = 0;
                two.top_left.x = 0;
                two.top_left.y = y_size;
                two.bottom_right = clust.areas[1].bottom_right;
                three.top_left.x = clust.areas[0].top_left.x;
                three.top_left.y = y_size;
                three.bottom_right.x = x_size;
                three.bottom_right.y = clust.areas[0].bottom_right.y;
                clust.areas.clear();
                clust.areas.push_back(zero);
                clust.areas.push_back(one);
                clust.areas.push_back(two);
                clust.areas.push_back(three);
            }
            else {                                                              //split horizontally
                zero.top_left = clust.areas[0].top_left;
                zero.bottom_right.x = x_size;
                zero.bottom_right.y = 0;
                one.top_left.x = 0;
                one.top_left.y = clust.areas[0].top_left.y;
                one.bottom_right.x = clust.areas[0].bottom_right.x;
                one.bottom_right.y = 0;
                two.top_left.x = 0;
                two.top_left.y = y_size;
                two.bottom_right.x = clust.areas[0].bottom_right.x;
                two.bottom_right.y = clust.areas[1].bottom_right.y;
                three.top_left.x = clust.areas[1].top_left.x;
                three.top_left.y = y_size;
                three.bottom_right.x = x_size;
                three.bottom_right.y = clust.areas[1].bottom_right.y;
                clust.areas.clear();
                clust.areas.push_back(zero);
                clust.areas.push_back(one);
                clust.areas.push_back(two);
                clust.areas.push_back(three);
            }
        }
        else if (clust.areas.size() == 4){
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right.x = x_size;
            zero.bottom_right.y = 0;
            one.top_left.x = 0;
            one.top_left.y = clust.areas[1].top_left.y;
            one.bottom_right.x = clust.areas[1].bottom_right.x;
            one.bottom_right.y = 0;
            two.top_left.x = 0;
            two.top_left.y = y_size;
            two.bottom_right.x = clust.areas[2].bottom_right.x;
            two.bottom_right.y = clust.areas[2].bottom_right.y;
            three.top_left.x = clust.areas[3].top_left.x;
            three.top_left.y = y_size;
            three.bottom_right.x = x_size;
            three.bottom_right.y = clust.areas[3].bottom_right.y;
            clust.areas.clear();
            clust.areas.push_back(zero);
            clust.areas.push_back(one);
            clust.areas.push_back(two);
            clust.areas.push_back(three);
        }
    }
    else {
        if (clust.areas.size() == 1){
            AABB zero(Point(0,0), Point(0,0));
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right = clust.areas[0].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
        }
        else if (clust.areas.size() == 2){
            AABB zero(Point(0,0), Point(0,0));
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right = clust.areas[1].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
        }
        else if (clust.areas.size() == 4){
            AABB zero(Point(0,0), Point(0,0));
            zero.top_left = clust.areas[0].top_left;
            zero.bottom_right = clust.areas[2].bottom_right;
            clust.areas.clear();
            clust.areas.push_back(zero);
        }
    }
}

AABB testFindSearchRange(AABB area, double step_len){
    AABB range = AABB(Point(0,0), Point(0,0));
    range.top_left.x = area.top_left.x - step_len;
    range.top_left.y = area.top_left.y + step_len;
    range.bottom_right.x = area.bottom_right.x + step_len;
    range.bottom_right.y = area.bottom_right.y - step_len;
    return range;
}

std::vector<Cluster> testBPcolCheck(AABB search_range, Quadtree &tree){
    std::vector<Cluster> ret_vec;
    ret_vec = tree.queryRange(search_range);
//    std::cout << "ret_vec.size() (BPColCheck) = " << ret_vec.size() << std::endl;
    if (ret_vec.size() == 1){
//        std::cout << "ret_vec[0].index (BPColCheck) = " << ret_vec[0].index << std::endl;
    }
    if (ret_vec.size() > 1){
        return ret_vec;
    }
    else {
        return {};
    }
}


int nbrAreasOut(Cluster* clust, Cluster* other, int x_size, int y_size){
    if ((clust != nullptr) && (other != nullptr)){
        if ((crossXBound(*clust, x_size)) && (crossYBound(*other, y_size))){
            return 4;
        }
        else if((crossYBound(*clust, y_size)) && (crossXBound(*other, x_size))){
            return 4;
        }
        else if (clust->areas.size() <= other->areas.size()){
            return other->areas.size();
        }
        else {
            return clust->areas.size();
        }
    }
    else {
        return 0;
    }
}

bool otherLeftTwo(Cluster* clust, Cluster* other){
    if (clust->areas[1].bottom_right.x - other->areas[0].top_left.x < 0.1){
        return true;
    }
    else {
        return false;
    }
}

bool clustLeftTwo(Cluster* clust, Cluster* other){
    if (clust->areas[0].top_left.x - other->areas[1].bottom_right.x < 0.1){
        return true;
    }
    else {
        return false;
    }
}

bool otherRightTwo(Cluster* clust, Cluster* other){
    if (clust->areas[0].top_left.x - other->areas[0].bottom_right.x < 0.1){
        return true;
    }
    else {
        return false;
    }
}

bool clustRightTwo(Cluster* clust, Cluster* other){
    if (other->areas[0].top_left.x - clust->areas[0].bottom_right.x < 0.1){
        return true;
    }
    else {
        return false;
    }
}

bool otherUpTwo(Cluster* clust, Cluster* other){
    if (other->areas[0].top_left.y - clust->areas[1].bottom_right.y < 0.1){
        return true;
    }
    else {
        return false;
    }
}

bool otherDownTwo(Cluster* clust, Cluster* other){
    if (clust->areas[0].top_left.y - other->areas[0].bottom_right.y < 0.1){
        return true;
    }
    else {
        return false;
    }
}

bool otherLeftFour(Cluster* clust, Cluster* other){
    if (clust->areas.size() == 4){
        if (other->areas.size() == 1){
            if (other->areas[0].top_left.x - clust->areas[1].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
        else if (other->areas.size() == 2){
            if (other->areas[0].top_left.x - clust->areas[1].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
    }
    else if (other->areas.size() == 4){
        if (clust->areas.size() == 1){
            if (clust->areas[0].top_left.x - other->areas[1].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
        else if (clust->areas.size() == 2){
            if (clust->areas[0].top_left.x - other->areas[1].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
    }
    return false;
}

bool otherRightFour(Cluster* clust, Cluster* other){
    if (clust->areas.size() == 4){
        if (other->areas.size() == 1){
            if (clust->areas[0].top_left.x - other->areas[0].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
        else if (other->areas.size() == 2){
            if (clust->areas[0].top_left.x - other->areas[0].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
    }
    else if (other->areas.size() == 4){
        if (clust->areas.size() == 1){
            if (other->areas[0].top_left.x - clust->areas[0].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
        else if (clust->areas.size() == 2){
            if (clust->areas[0].top_left.x - other->areas[0].bottom_right.x < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
    }
    return false;
}

bool otherDownFour(Cluster* clust, Cluster* other){
    if (clust->areas.size() == 4){
        if (other->areas.size() == 1){
            if (other->areas[0].bottom_right.y - clust->areas[1].top_left.y < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
        else if (other->areas.size() == 2){
            //Write this when done with the simple cases.
        }
    }
    else if (other->areas.size() == 4){
        if (clust->areas.size() == 1){
            if (clust->areas[0].bottom_right.y - other->areas[1].top_left.y < 0.1){
                return true;
            }
            else {
                return false;
            }
            //Write this when done with the simple cases.
        }
        else if (clust->areas.size() == 2){

        }
    }
    return false;

}

bool otherUpFour(Cluster* clust, Cluster* other){
    if (clust->areas.size() == 4){
        if (other->areas.size() == 1){
            if (clust->areas[2].bottom_right.y - other->areas[0].top_left.y < 0.1){
                return true;
            }
            else {
                return false;
            }
        }
        else if (other->areas.size() == 2){
            //Write this when done with the simple cases.
        }
    }
    else if (other->areas.size() == 4){
        if (clust->areas.size() == 1){
            if (other->areas[2].bottom_right.y - clust->areas[0].top_left.y < 0.1){
                return true;
            }
            else {
                return false;
            }
            //Write this when done with the simple cases.
        }
        else if (clust->areas.size() == 2){

        }
    }
    return false;

}


void updateAreas(Cluster* clust, Cluster* other, int x_size, int y_size){
    if (nbrAreasOut(clust, other, x_size, y_size) == 1){
        clust->areas[0].top_left.x = std::min(clust->areas[0].top_left.x,
                                              other->areas[0].top_left.x);
        clust->areas[0].top_left.y = std::max(clust->areas[0].top_left.y,
                                              other->areas[0].top_left.y);
        clust->areas[0].bottom_right.x=std::max(clust->areas[0].bottom_right.x,
                                                other->areas[0].bottom_right.x);
        clust->areas[0].bottom_right.y=std::min(clust->areas[0].bottom_right.y,
                                                other->areas[0].bottom_right.y);
        return;
    }
    else if (nbrAreasOut(clust, other, x_size, y_size) == 2){
        if (crossXBound(*clust, x_size)){
            if (crossXBound(*other, x_size)){                                    //This checks 3) from the notes, where both rectangles crosses the boundary
                clust->areas[0].top_left.x = std::min(clust->areas[0].top_left.x,
                                                      other->areas[0].top_left.x);
                clust->areas[0].top_left.y = std::max(clust->areas[0].top_left.y,
                                                      other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x = x_size;
                clust->areas[0].bottom_right.y = std::min(clust->areas[0].bottom_right.y,
                                                          other->areas[0].bottom_right.y);

                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y = std::max(clust->areas[1].top_left.y,
                                                      other->areas[1].top_left.y);
                clust->areas[1].bottom_right.x=std::max(clust->areas[1].bottom_right.x,
                                                        other->areas[1].bottom_right.x);
                clust->areas[1].bottom_right.y = clust->areas[0].bottom_right.y;
                return;
            }
            else if (otherRightTwo(clust, other)){                              //This checks 1) from the notes
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x = x_size;
                clust->areas[0].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[0].bottom_right.y);
                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y=std::max(clust->areas[1].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[1].bottom_right.x = clust->areas[1].bottom_right.x;
                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[0].bottom_right.y);
                return;
            }
            else if (otherLeftTwo(clust, other)){                               //This checks 2) from the notes
                clust->areas[0].top_left.x = clust->areas[0].top_left.x;
                clust->areas[0].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x = x_size;
                clust->areas[0].bottom_right.y =
                        std::min(clust->areas[0].bottom_right.y,
                                 other->areas[0].bottom_right.y);
                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y=std::max(clust->areas[1].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[1].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[0].bottom_right.y);
                return;
            }
            else {
                std::cout << "There is an error in updateAreas (A)" << std::endl;
            }
        }
        else if (crossXBound(*other, x_size)){                                  //The other cluster crosses the boundary. Have to expand the Areas vector.
            clust->areas.push_back(clust->areas[0]);                            //This is done by just duplicating the 0th element.
            if (clustRightTwo(clust, other)){                                   //This checks 4) in the notes.
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x = x_size;
                clust->areas[0].bottom_right.y =
                        std::min(clust->areas[0].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y=std::min(clust->areas[0].top_left.y,
                                                    other->areas[1].top_left.y);
                clust->areas[1].bottom_right.x = other->areas[1].bottom_right.x;
                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[0].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                return;
            }
            else if (clustLeftTwo(clust, other)){                               //checks 5) in the notes.
                clust->areas[0].top_left.x = other->areas[0].top_left.x;
                clust->areas[0].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x = x_size;
                clust->areas[0].bottom_right.y =
                        std::min(clust->areas[0].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[1].top_left.y);
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[1].bottom_right.x,
                                 other->areas[1].bottom_right.x);
                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                return;
            }
            else {
                std::cout << "There is an error in updateAreas (B)" << std::endl;
            }
        }
        else if(crossYBound(*clust, y_size)){
            if (crossYBound(*other, y_size)){                                   //checks 6) in the notes
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[0].bottom_right.y = 0;
                clust->areas[1].top_left.x=std::min(clust->areas[1].top_left.x,
                                                    other->areas[1].top_left.x);
                clust->areas[1].top_left.y = y_size;
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[1].bottom_right.x,
                                 other->areas[1].bottom_right.x);
                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                return;
            }
            else if (otherUpTwo(clust, other)){                                 //checks 7) in the notes
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y =std::max(clust->areas[0].top_left.y,
                                                     other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[0].bottom_right.y = 0;
                clust->areas[1].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[1].top_left.y = y_size;
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[1].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[1].bottom_right.y =
                        std::max(clust->areas[1].bottom_right.y,
                                 other->areas[0].bottom_right.y);
                return;
            }
            else if (otherDownTwo(clust, other)){                               //checks 8) in the notes.
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y=std::min(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[0].bottom_right.y = 0;
                clust->areas[1].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[1].top_left.y = y_size;
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[1].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[0].bottom_right.y);
                return;
            }
            else {
                std::cout << "There is an error in updateAreas (C)" << std::endl;
            }
        }
        else if (crossYBound(*other, y_size)){
            clust->areas.push_back(clust->areas[0]);                            //We have to expand the areas vector.
            if (otherUpTwo(other, clust)){                                      //checks 9) in notes. Note swapped argument
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y=std::max(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[0].bottom_right.y = 0;
                clust->areas[1].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[1].top_left.y = y_size;
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[1].bottom_right.x);
                clust->areas[1].bottom_right.y =
                        std::max(clust->areas[1].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                return;
            }
            else if (otherDownTwo(other, clust)){                               //checks 10) in the notes. Note swapped argument
                clust->areas[0].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[0].top_left.x);
                clust->areas[0].top_left.y=std::min(clust->areas[0].top_left.y,
                                                    other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[0].bottom_right.x);
                clust->areas[0].bottom_right.y = 0;
                clust->areas[1].top_left.x=std::min(clust->areas[0].top_left.x,
                                                    other->areas[1].top_left.x);
                clust->areas[1].top_left.y = y_size;
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[0].bottom_right.x,
                                 other->areas[1].bottom_right.x);

                clust->areas[1].bottom_right.y =
                        std::min(clust->areas[1].bottom_right.y,
                                 other->areas[1].bottom_right.y);
                return;
            }
            else {
                std::cout << "There is an error in updateAreas (D)" << std::endl;
            }
        }
    }
    else if (nbrAreasOut(clust, other, x_size, y_size) == 4){
        if ((crossXBound(*clust, x_size)) && (crossYBound(*clust, y_size)) &&
             (other->areas.size() == 1)){
            if (otherLeftFour(clust, other)){
                if (otherDownFour(clust, other)){                               //checks 11) in the notes
                    clust->areas[0].top_left.y =
                            std::max(clust->areas[1].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[1].top_left.x = 0;
                    clust->areas[1].top_left.y =
                            std::max(clust->areas[1].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[1].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[0].bottom_right.x);
                    clust->areas[2].top_left.x = 0;
                    clust->areas[2].top_left.y = y_size;
                    clust->areas[2].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[0].bottom_right.x);
                    clust->areas[2].bottom_right.y =
                            clust->areas[2].bottom_right.y;
                    return;

                }
                else if (otherUpFour(clust, other)){                            //checks 12) in notes
                    clust->areas[1].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[0].bottom_right.x);
                    clust->areas[2].bottom_right.x =
                            std::max(clust->areas[2].bottom_right.x,
                                     other->areas[0].bottom_right.x);
                    clust->areas[2].bottom_right.y =
                            std::min(clust->areas[2].bottom_right.y,
                                     other->areas[0].bottom_right.y);
                    clust->areas[3].bottom_right.y =
                            std::min(clust->areas[3].bottom_right.y,
                                     other->areas[0].bottom_right.y);
                    return;
                }
                else {
                    std::cout << "There is an error in updateAreas (E)" << std::endl;
                }
            }
            else if (otherRightFour(clust,other)){
                if (otherUpFour(clust, other)){                               //checks 13) in the notes
                    clust->areas[0].top_left.x =
                            std::min(clust->areas[0].top_left.x,
                                     other->areas[0].top_left.x);
                    clust->areas[2].bottom_right.y =
                            std::min(clust->areas[2].bottom_right.y,
                                     other->areas[0].bottom_right.y);
                    clust->areas[3].top_left.x =
                            std::min(clust->areas[3].top_left.x,
                                     other->areas[0].top_left.x);
                    clust->areas[3].bottom_right.y =
                            std::min(clust->areas[3].bottom_right.y,
                                     other->areas[0].bottom_right.y);
                    return;
                }
                else if (otherDownFour(clust, other)){                          //checks 14) in the notes
                    clust->areas[0].top_left.x =
                            std::min(clust->areas[0].top_left.x,
                                     other->areas[0].top_left.x);
                    clust->areas[0].top_left.y =
                            std::max(clust->areas[0].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[1].top_left.y =
                            std::max(clust->areas[1].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[3].top_left.x =
                            std::min(clust->areas[3].top_left.x,
                                     other->areas[0].top_left.x);
                    return;
                }
                else {
                    std::cout << "There is an error in updateAreas (F)" << std::endl;
                }
            }
        }
        else if ((crossXBound(*other, x_size)) && (crossYBound(*other, x_size))
                 && (clust->areas.size() == 1)){
            if (otherLeftFour(clust, other)){
                if (otherDownFour(clust, other)){                               //checks 15) in the notes
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas[0].top_left.x = other->areas[0].top_left.x;
                    clust->areas[0].top_left.y =
                            std::max(other->areas[0].top_left.y,
                                     clust->areas[0].top_left.y);
                    clust->areas[0].bottom_right = other->areas[0].bottom_right;
                    clust->areas[1].top_left.x = 0;
                    clust->areas[1].top_left.y =
                            std::max(other->areas[1].top_left.y,
                                     clust->areas[1].top_left.y);
                    clust->areas[1].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[1].bottom_right.x);
                    clust->areas[1].bottom_right.y = 0;
                    clust->areas[2].top_left = other->areas[2].top_left;
                    clust->areas[2].bottom_right.x =
                            std::max(clust->areas[2].bottom_right.x,
                                     other->areas[2].bottom_right.x);
                    clust->areas[2].bottom_right.y =
                            other->areas[2].bottom_right.y;
                    clust->areas[3].top_left = other->areas[3].top_left;
                    clust->areas[3].bottom_right = other->areas[3].bottom_right;
                    return;
                }
                else if (otherUpFour(clust, other)){                            //checks 16) in the notes
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas[0].top_left = other->areas[0].top_left;
                    clust->areas[0].bottom_right = other->areas[0].bottom_right;
                    clust->areas[1].top_left = other->areas[1].top_left;
                    clust->areas[1].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[1].bottom_right.x);
                    clust->areas[1].bottom_right.y = 0;
                    clust->areas[2].top_left = other->areas[2].top_left;
                    clust->areas[2].bottom_right.x =
                            std::max(clust->areas[2].bottom_right.x,
                                     other->areas[2].bottom_right.x);
                    clust->areas[2].bottom_right.y =
                            std::min(clust->areas[2].bottom_right.y,
                                     other->areas[2].bottom_right.y);
                    clust->areas[3].top_left = other->areas[3].top_left;
                    clust->areas[3].bottom_right.x = x_size;
                    clust->areas[3].bottom_right.y =
                            std::min(clust->areas[3].bottom_right.y,
                                     other->areas[3].bottom_right.y);
                    return;
                }
                else {
                    std::cout << "There is an error in updateAreas (G)" << std::endl;
                }
            }
            else if (otherRightFour(clust, other)){
                if (otherUpFour(clust, other)){                                 //checks 17 in the notes
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas[0].top_left.x =
                            std::min(clust->areas[0].top_left.x,
                                     other->areas[0].top_left.x);
                    clust->areas[0].top_left.y = other->areas[0].top_left.y;
                    clust->areas[0].bottom_right = other->areas[0].bottom_right;
                    clust->areas[1].top_left = other->areas[1].top_left;
                    clust->areas[1].bottom_right = other->areas[1].bottom_right;
                    clust->areas[2].top_left = other->areas[2].top_left;
                    clust->areas[2].bottom_right.x =
                            other->areas[2].bottom_right.x;
                    clust->areas[2].bottom_right.y =
                            std::min(clust->areas[2].bottom_right.y,
                                     other->areas[2].bottom_right.y);
                    clust->areas[3].top_left.x =
                            std::min(clust->areas[3].top_left.x,
                                     other->areas[3].top_left.x);
                    clust->areas[3].top_left.y = other->areas[3].top_left.y;
                    clust->areas[3].bottom_right.x = x_size;
                    clust->areas[3].bottom_right.y =
                            std::min(clust->areas[3].bottom_right.y,
                                     other->areas[3].bottom_right.y);
                    return;
                }
                else if (otherDownFour(clust, other)){                          //checks 18 in the notes
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas.push_back(clust->areas[0]);
                    clust->areas[0].top_left.y =
                            std::max(clust->areas[0].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[0].top_left.x =
                            std::min(clust->areas[0].top_left.x,
                                     other->areas[0].top_left.x);
                    clust->areas[0].bottom_right = other->areas[0].bottom_right;
                    clust->areas[1].top_left.x = other->areas[0].top_left.x;
                    clust->areas[1].top_left.y =
                            std::max(clust->areas[1].top_left.y,
                                     other->areas[1].top_left.y);
                    clust->areas[1].bottom_right = other->areas[1].bottom_right;
                    clust->areas[2] = other->areas[2];
                    clust->areas[3].top_left.x =
                            std::min(clust->areas[3].top_left.x,
                                     other->areas[3].top_left.x);
                    clust->areas[3].top_left.y = y_size;
                    clust->areas[3].bottom_right = other->areas[3].bottom_right;
                    return;
                }
                else {
                    std::cout << "There is an error in updateAreas (H)" << std::endl;
                }
            }
        }
        else if ((crossXBound(*clust, x_size)) && (crossYBound(*clust, y_size))
                  && (other->areas.size() == 2)){
            if (crossYBound(*other, y_size)){
                if (otherRightFour(clust, other)){                               //checks 19) in the notes
                    clust->areas[0].top_left.x =
                            std::min(clust->areas[0].top_left.x,
                                     other->areas[0].top_left.x);
                    clust->areas[0].top_left.y =
                            std::max(clust->areas[0].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[3].top_left.x =
                            std::min(clust->areas[3].top_left.x,
                                     other->areas[1].top_left.x);
                    clust->areas[3].bottom_right.y =
                            std::min(clust->areas[3].bottom_right.y,
                                     other->areas[1].bottom_right.y);
                    return;
                }
                else if (otherLeftFour(clust, other)){                          //checks 20) in the notes
                    clust->areas[1].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[0].bottom_right.x);
                    clust->areas[1].top_left.y =
                            std::max(clust->areas[1].top_left.y,
                                     other->areas[0].top_left.y);
                    clust->areas[2].bottom_right.x =
                            std::max(clust->areas[1].bottom_right.x,
                                     other->areas[0].bottom_right.x);
                    clust->areas[2].bottom_right.y =
                            std::min(clust->areas[2].bottom_right.y,
                                     other->areas[1].bottom_right.y);
                    return;
                }
                else {
                    std::cout << "There is an error in updateAreas (I)" << std::endl;
                }
            }                                                                   //case 21) and 22) are mapped to 11) and 12), but works.
        }
        else if ((crossXBound(*other, x_size)) && (crossYBound(*other, y_size))
                 && (clust->areas.size() == 2)){
            if (otherLeftFour(clust, other)){                                   //checks 23) in the notes
                clust->areas.push_back(clust->areas[0]);
                clust->areas.push_back(clust->areas[1]);
                clust->areas[0] = other->areas[0];
                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y =
                        std::min(clust->areas[0].top_left.y,
                                 other->areas[1].top_left.y);
                clust->areas[1].bottom_right.x =
                        std::max(clust->areas[2].bottom_right.x,
                                 other->areas[1].bottom_right.x);
                clust->areas[1].bottom_right.y = 0;
                clust->areas[2].top_left = other->areas[2].top_left;
                clust->areas[2].bottom_right.x =
                        std::max(clust->areas[3].bottom_right.x,
                                 other->areas[2].bottom_right.x);
                clust->areas[2].bottom_right.y =
                        std::min(clust->areas[3].bottom_right.y,
                                 other->areas[2].bottom_right.y);
                clust->areas[3] = other->areas[3];
                return;
            }
            else if (otherRightFour(clust, other)){                             //checks 24) in the notes
                clust->areas.push_back(clust->areas[0]);
                clust->areas.push_back(clust->areas[1]);
                clust->areas[0].top_left.x =
                        std::min(clust->areas[0].top_left.x,
                                 other->areas[0].top_left.x);
                clust->areas[0].top_left.y =
                        std::max(clust->areas[0].top_left.y,
                                 other->areas[0].top_left.y);
                clust->areas[0].bottom_right.x = x_size;
                clust->areas[0].bottom_right.y = 0;
                clust->areas[1].top_left.x = 0;
                clust->areas[1].top_left.y =
                        std::max(clust->areas[2].top_left.y,
                                 other->areas[1].top_left.y);
                clust->areas[1].bottom_right = other->areas[1].bottom_right;
                clust->areas[2].top_left = other->areas[2].top_left;
                clust->areas[2].bottom_right.x = other->areas[2].bottom_right.x;
                clust->areas[2].bottom_right.y =
                        std::min(clust->areas[3].bottom_right.y,
                                 other->areas[2].bottom_right.y);
                clust->areas[3].bottom_right.x = x_size;
                clust->areas[3].bottom_right.y =
                        std::min(clust->areas[3].bottom_right.y,
                                 other->areas[3].bottom_right.y);
                return;
            }
            else {
                std::cout << "There is an error in updateAreas (J)" << std::endl;
            }
        }
    }
}

double calcStepLen(Cluster* clust, double gamma, double PI, double D_0,
                   double dt){
    if (clust != nullptr){
        double sum = 0;
        int counter = 0;
        for (auto&& particle:clust->particles){
            sum += PI*particle.r_p*particle.r_p;
            counter ++;
        }
        return std::sqrt(2*dt*D_0*std::pow(sum, -gamma));
    }
}








bool consMomentum(int clustIdx, int otherIdx){
    if (clustIdx < otherIdx){
        return true;
    }
    else {
        return false;
    }
}






