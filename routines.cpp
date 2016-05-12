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
    int org_size = clust_test.size();
    int part_placed = clust_test.size();
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;
    velocity vel;

    while (int(clust_test.size()) < org_size + nbr_particles){
        tmp.x = coord_x();
        tmp.y = coord_y();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clust_test.size()); ++i){
            dist = findDistance(clust_test[i]->particles[0].pos, tmp, len, height);
            if (clust_test[i]->particles[0].r_p + r_p > dist){
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
            clust_test.insert(std::make_pair(part_placed,
                                             new Cluster(false, areas_tmp,
                                                         tmp_part, part_placed,
                                                         r_p, vel, tmp)));
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;                                                       //We reset and prepare to place the next particle.
    }
}

void distributeDust(int len, int height, int nbr_particles,
                    std::mt19937::result_type seed_x, bool varying_size,
                    std::mt19937::result_type seed_y, double r_p,
                    std::map<int, Cluster*> &clusters){
    auto coord_x = std::bind(std::uniform_real_distribution<double>(0,len),     //We define a function to generate a random point in the x-dimension on the domain
                               std::mt19937(seed_x));
    auto coord_y = std::bind(std::uniform_real_distribution<double>(0,height),  //similarly for the y-dimension.
                               std::mt19937(seed_y));
    auto rand_size = std::bind(std::uniform_real_distribution<float>(1.0,5.0),
                               std::mt19937(seed_x));
    double dist;
    int container_position = clusters.size();
    int org_size = container_position;
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;
    velocity vel;

    while (int(clusters.size()) < org_size + nbr_particles){
        tmp.x = coord_x();
        tmp.y = coord_y();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clusters.size()); ++i){
            dist = findDistance(clusters[i]->particles[0].pos, tmp, len, height);
            if (clusters[i]->particles[0].r_p + r_p > dist){
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
            clusters.insert(std::make_pair(container_position,
                                             new Cluster(false, areas_tmp,
                                                         tmp_part,
                                                         container_position,
                                                         r_p, vel, tmp)));
            std::cout << "dust at: " << container_position << std::endl;
            container_position++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;                                                       //We reset and prepare to place the next particle.
    }
}

void distributeDustTestSection(int len, int height, int nbr_particles,
                    std::mt19937::result_type seed_x, bool varying_size,
                    std::mt19937::result_type seed_y, double r_p,
                    std::map<int, Cluster*> &clusters){
    auto coord_x = std::bind(std::uniform_real_distribution<double>(len/2.0 - len/4.0, len/2.0 + len/4.0),     //We define a function to generate a random point in the x-dimension on the domain
                               std::mt19937(seed_x));
    auto coord_y = std::bind(std::uniform_real_distribution<double>(height/2.0 - height/4.0, height/2.0 + height/4.0),  //similarly for the y-dimension.
                               std::mt19937(seed_y));
    auto rand_size = std::bind(std::uniform_real_distribution<float>(1.0,5.0),
                               std::mt19937(seed_x));
    double dist;
    int container_position = clusters.size();
    int org_size = container_position;
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;
    velocity vel;

    while (int(clusters.size()) < org_size + nbr_particles){
        tmp.x = coord_x();
        tmp.y = coord_y();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clusters.size()); ++i){
            dist = findDistance(clusters[i]->particles[0].pos, tmp, len, height);
            if (clusters[i]->particles[0].r_p + r_p > dist){
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
            clusters.insert(std::make_pair(container_position,
                                             new Cluster(false, areas_tmp,
                                                         tmp_part,
                                                         container_position,
                                                         r_p, vel, tmp)));
            std::cout << "dust at: " << container_position << std::endl;
            container_position++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;                                                       //We reset and prepare to place the next particle.
    }
}

void distributeParticlesTestSection(int len, int height, int nbr_particles,
                            std::vector<Cluster> &clusters,
                            std::mt19937::result_type seed_x,
                            std::mt19937::result_type seed_y, double r_p,
                            bool varying_size, std::map<int, Cluster*> &clust_test){
    auto coord_x = std::bind(std::uniform_real_distribution<double>(len/2.0 - len/4.0, len/2.0 + len/4.0),     //We define a function to generate a random point in the x-dimension on the domain
                               std::mt19937(seed_x));
    auto coord_y = std::bind(std::uniform_real_distribution<double>(height/2.0 - height/4.0, height/2.0 + height/4.0),  //similarly for the y-dimension.
                               std::mt19937(seed_y));
    auto rand_size = std::bind(std::uniform_real_distribution<float>(1.0,5.0),
                               std::mt19937(seed_x));
    double dist;
    int org_size = clust_test.size();
    int part_placed = clust_test.size();
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;
    velocity vel;

    while (int(clust_test.size()) < org_size + nbr_particles){
        tmp.x = coord_x();
        tmp.y = coord_y();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clust_test.size()); ++i){
            dist = findDistance(clust_test[i]->particles[0].pos, tmp, len, height);
            if (clust_test[i]->particles[0].r_p + r_p > dist){
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
            clust_test.insert(std::make_pair(part_placed,
                                             new Cluster(false, areas_tmp,
                                                         tmp_part, part_placed,
                                                         r_p, vel, tmp)));
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

void testPutInQuadtree(std::map<int, Cluster*> clusters, Quadtree &qtree){
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end(); iterator++){
        if (iterator->second != nullptr){
            qtree.insert(*(iterator->second));
        }
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

std::vector<Cluster> BPcolCheck(Cluster* cluster,
                                std::vector<AABB> search_ranges,
                                Quadtree &tree){
    std::vector<Cluster> ret_vec;
//    if (cluster != nullptr){
        std::vector<Cluster> tmp_res;
        int index = cluster->index;
        for (auto&& search_range:search_ranges){
            tmp_res = tree.queryRange(search_range);
//            if (tmp_res.size() > 1){
                for (auto&& res: tmp_res){
//                    if (res.index != index){
                        ret_vec.push_back(res);
//                    }
//                }
//            }
        }
    }
    return ret_vec;
}

double LHit(double step_L, double step_dir, Particle one, Particle two,
            int x_size, int y_size){
    //rewrite the check to take findDistance into account!
    step_L += 0.001;
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

//    double distance = std::sqrt(dx*dx + dy*dy) - one.r_p - two.r_p;
//    if (distance < step_L){
//        return distance;
//    }
//    else {
//        return step_L;
//    }


    double d_p = one.r_p + two.r_p;
    double a = 1.0;
    double b = 2.0*(std::cos(step_dir)*(dx) +
                    std::sin(step_dir)*(dy));
    double c = (-dx)*(-dx) +
               (-dy)*(-dy) - d_p*d_p;
    double res[2];
    bool sol = true;
    if ((b*b - 4*a*c) < 0){
        sol = false;
        return step_L;
    }
    else {
        res[0] = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
        res[1] = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
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

double NPColCheckOrg(Cluster* cluster, std::vector<Cluster> targets,
                  double step_len, double step_dir, int &col_with,
                  int x_size, int y_size){
    double tmp;
    double ret_val = step_len;
    if (cluster != nullptr){
        for (auto&& target:targets){
            if (target.index != cluster->index){
                for (auto&& particle:cluster->particles){
                    for (auto&& tar_part: target.particles){
                        tmp = LHit(step_len, step_dir, particle, tar_part,
                                   x_size, y_size);
                        if (tmp < ret_val){
                            ret_val = tmp;
                            col_with = target.index;
                        }
                    }
                }
            }
        }
    }
    return ret_val;                                                             //Returns the correct step length
}


double NPColCheck(Cluster* cluster, std::vector<Cluster> targets,
                  double &step_len, double step_dir, int &col_with,
                  int x_size, int y_size, Quadtree_vel &vel_tree, double dt,
                  double L_typical, double rho_air, double rho_dust,
                  double C_sphere, double PI, bool &will_collide){
    double tmp;
    double ret_val;
    double sx = 0, sy = 0;
    if (cluster != nullptr){
        will_collide = false;
        ret_val = step_len;
        for (auto&& target : targets){
            if (target.index != cluster->index){
                for (auto&& particle: cluster->particles){
                    for (auto&& tar_part: target.particles){
                        tmp = LHit(step_len, step_dir, particle, tar_part,
                                   x_size, y_size);
                        if (std::abs(tmp - ret_val) < 0.00001){
                            will_collide = true;
                            col_with = target.index;
                        }
                        if (tmp < ret_val){
                            ret_val = tmp;
                            col_with = target.index;
                        }
                    }
                }
            }
        }
    }
    return ret_val;                                                             //Returns the correct step length
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
                      int y_size, double rho_carbon, double rho_dust,
                      double r_dust, double r_carbon, double PI,
                      double L_typical){
    if ((clust != nullptr) && (other != nullptr)){
        if (clust->index != other->index){
            for (auto&& particle:clust->particles){
                if (particle.r_p > 1000){
                    std::cout << "something is wrong cluster" << std::endl;
                    std::cout << particle.r_p << std::endl;
                    std::cout << "particle.size() = " << clust->particles.size() << std::endl;

                }
            }
            for (auto&& particle:other->particles){
                if (particle.r_p > 1000){
                    particle.r_p = 5.0;
                }
            }
            int org_size = clust->particles.size();
            updateAreas(clust, other, x_size, y_size);
            clust->particles.insert(clust->particles.end(),
                                    other->particles.begin(),
                                    other->particles.end());
            if (clust->vel.v > other->vel.v){
                clust->vel = other->vel;
            }
            clust->CM = FindCM(*clust, rho_carbon, rho_dust, r_dust, r_carbon,
                   PI, L_typical);
            double diameter = 0;
            for (int i = 0; i < int(clust->particles.size()); ++i){
                for (int j = i+1; j < int(clust->particles.size()); ++j){
                    double tmp_diameter = findDistance(clust->particles[i].pos,
                                                       clust->particles[j].pos,
                                                       x_size, y_size) +
                                          clust->particles[i].r_p +
                                          clust->particles[j].r_p;
                    if (tmp_diameter> diameter){
                        diameter = tmp_diameter;
                    }
                }
            }
            clust->radius = diameter/2.0;
            if (clust->particles.size() != org_size + other->particles.size()){
                std::cout << "size missmatch for particles" << std::endl;
            }
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
    if (clust->areas[1].bottom_right.x > other->areas[0].top_left.x) {
        return true;
    }
    else {
        return false;
    }
}

bool testOtherLeft(Cluster* other, int x_size){
    if (other->areas.size() == 1){
        if ((other->areas[0].top_left.x < x_size/2.0) ||
            (other->areas[0].bottom_right.x < x_size/2.0)){
//            std::cout << "other->areas[0].top_left.x = " << (other->areas[0].top_left.x) << std::endl;
//            std::cout << " x_size/2.0 = " <<  x_size/2.0 << std::endl;
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

bool testOtherRight(Cluster* other, int x_size){
    if (other->areas.size() == 1){
        if ((other->areas[0].top_left.x > x_size/2.0) ||
            (other->areas[0].bottom_right.x > x_size/2.0)){
//            std::cout << "other->areas[0].top_left.x = " << (other->areas[0].top_left.x) << std::endl;
//            std::cout << " x_size/2.0 = " <<  x_size/2.0 << std::endl;
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

bool testClustLeft(Cluster* clust, int x_size){
//    std::cout << "clust->areas[0].top_left.x = " << (clust->areas[0].top_left.x) << std::endl;
//    std::cout << "clust->areas[0].bottom_right.x = " << (clust->areas[0].bottom_right.x) << std::endl;
//    std::cout << " x_size/2.0 = " <<  x_size/2.0 << std::endl;
    if (clust->areas.size() == 1){
        if ((clust->areas[0].top_left.x < x_size/2.0) ||
            (clust->areas[0].bottom_right.x < x_size/2.0)){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust->areas.size() == 2){
        if ((clust->areas[0].top_left.x < x_size/2.0) ||
            (clust->areas[0].bottom_right.x < x_size/2.0)){
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

bool testClustRight(Cluster* clust, int x_size){
    if (clust->areas.size() == 1){
        if (clust->areas[0].top_left.x > x_size/2.0){
            return true;
        }
        else if (clust->areas[0].bottom_right.x > x_size/2.0){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust->areas.size() == 2){
        if (clust->areas[0].top_left.x > x_size/2.0){
            return true;
        }
        else if (clust->areas[0].bottom_right.x > x_size/2.0){
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

bool testClustUp(Cluster* clust, int y_size){
    if (clust->areas.size() == 1){
        if (clust->areas[0].top_left.y > y_size/2.0){
            return true;
        }
        else if (clust->areas[0].bottom_right.y > y_size/2.0){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust->areas.size() == 2){
        if (clust->areas[0].top_left.y > y_size/2.0){
            return true;
        }
        else if (clust->areas[0].bottom_right.y > y_size/2.0){
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

bool testClustDown(Cluster* clust, int y_size){
    if (clust->areas.size() == 1){
        if (clust->areas[0].top_left.y < y_size/2.0){
            return true;
        }
        else if (clust->areas[0].bottom_right.y < y_size/2.0){
            return true;
        }
        else {
            return false;
        }
    }
    else if (clust->areas.size() == 2) {
        if (clust->areas[0].top_left.y < y_size/2.0){
            return true;
        }
        else if (clust->areas[0].bottom_right.y < y_size/2.0){
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

bool clustLeftTwo(Cluster* clust, Cluster* other){
    if (clust->areas[0].top_left.x < other->areas[1].bottom_right.x){
        return true;
    }
    else {
        return false;
    }
}

bool otherRightTwo(Cluster* clust, Cluster* other){
    if (clust->areas[1].top_left.x > other->areas[0].top_left.x){
        return true;
    }
    else {
        return false;
    }
}

bool clustRightTwo(Cluster* clust, Cluster* other){
    if (other->areas[0].top_left.x < clust->areas[0].bottom_right.x){
        return true;
    }
    else {
        return false;
    }
}

bool otherUpTwo(Cluster* clust, Cluster* other){
    if (other->areas[0].top_left.y > clust->areas[1].bottom_right.y){
        return true;
    }
    else {
        return false;
    }
}

bool otherDownTwo(Cluster* clust, Cluster* other){
    if (clust->areas[0].top_left.y < other->areas[0].bottom_right.y){
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
    if (nbrAreasOut(clust, other, x_size, y_size) == 3){
        std::cout << "hi" << std::endl;
        std::cout << "clust areas size: " << clust->areas.size() << std::endl;
        std::cout << "other areas size: " << other->areas.size() << std::endl;
    }
//    std::cout << clust->index << " has collided with " << other->index << std::endl;
//    std::cout << "clust coord: " << clust->particles[0].pos.x << "," << clust->particles[0].pos.y << std::endl;
//    std::cout << "other coord: " << other->particles[0].pos.x << "," << other->particles[0].pos.y << std::endl;

//    std::cout << "clust[0] tlx = " << clust->areas[0].top_left.x << std::endl;
//    std::cout << "clust[0] brx = " << clust->areas[0].bottom_right.x << std::endl;
//    std::cout << "other[0] tlx = " << other->areas[0].top_left.x << std::endl;
//    std::cout << "other[0] brx = " << other->areas[0].bottom_right.x << std::endl;


//    std::cout << "nbr areas out = " << nbrAreasOut(clust, other, x_size, y_size) << std::endl;
//    std::cout << "does the cluster we focus on cross the x-boundary? " << crossXBound(*clust, x_size) << std::endl;
//    std::cout << "does the target cross the x-boundary? " << crossXBound(*other, x_size) << std::endl;
//    std::cout << "does the cluster we focus on cross the y-boundary? " << crossYBound(*clust, x_size) << std::endl;
//    std::cout << "does the target cross the y-boundary? " << crossYBound(*other, x_size) << std::endl;
//    if (nbrAreasOut(clust, other, x_size, y_size) == 2){
//        std::cout << "is the target on the left side of the domain? " << testOtherLeft(other, x_size) << std::endl;
//        std::cout << "is the target on the right side of the domain? " << testOtherRight(other, x_size) << std::endl;
//        std::cout << "is the cluster on the left side of the domain? " << testClustLeft(clust, x_size) << std::endl;
//        std::cout << "is the cluster on the right side of the domain? " << testClustRight(clust, x_size) << std::endl;
//        std::cout << "is the target on the top side of the domain? " << testClustUp(other, y_size) << std::endl;
//        std::cout << "is the cluster on the top side of the domain? " <<testClustUp(clust, y_size) << std::endl;
//    }
//    std::cout << std::endl;


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
            else if (testOtherRight(other, x_size)){                              //This checks 1) from the notes
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
            else if (testOtherLeft(other, x_size)){                               //This checks 2) from the notes
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
                std::cout << clust->index << " has collided with " << other->index << std::endl;
                std::cout << "clust coord: " << clust->particles[0].pos.x << "," << clust->particles[0].pos.y << std::endl;
                std::cout << "other coord: " << other->particles[0].pos.x << "," << other->particles[0].pos.y << std::endl;

                std::cout << "clust[0] tlx = " << clust->areas[0].top_left.x << std::endl;
                std::cout << "clust[0] brx = " << clust->areas[0].bottom_right.x << std::endl;
                std::cout << "other[0] tlx = " << other->areas[0].top_left.x << std::endl;
                std::cout << "other[0] brx = " << other->areas[0].bottom_right.x << std::endl;


                std::cout << "nbr areas out = " << nbrAreasOut(clust, other, x_size, y_size) << std::endl;
                std::cout << "does the cluster we focus on cross the x-boundary? " << crossXBound(*clust, x_size) << std::endl;
                std::cout << "does the target cross the x-boundary? " << crossXBound(*other, x_size) << std::endl;
                std::cout << "does the cluster we focus on cross the y-boundary? " << crossYBound(*clust, x_size) << std::endl;
                std::cout << "does the target cross the y-boundary? " << crossYBound(*other, x_size) << std::endl;
                if (nbrAreasOut(clust, other, x_size, y_size) == 2){
                    std::cout << "is the target on the left side of the domain? " << testOtherLeft(other, x_size) << std::endl;
                    std::cout << "is the target on the right side of the domain? " << testOtherRight(other, x_size) << std::endl;
                    std::cout << "is the cluster on the left side of the domain? " << testClustLeft(clust, x_size) << std::endl;
                    std::cout << "is the cluster on the right side of the domain? " << testClustRight(clust, x_size) << std::endl;
                    std::cout << "is the target on the top side of the domain? " << testClustUp(other, y_size) << std::endl;
                    std::cout << "is the cluster on the top side of the domain? " <<testClustUp(clust, y_size) << std::endl;
                }
                std::cout << std::endl;
            }
        }
        else if (crossXBound(*other, x_size)){                                  //The other cluster crosses the boundary. Have to expand the Areas vector.
            if (testClustRight(clust, x_size)){                                 //This checks 4) in the notes.
                clust->areas.push_back(clust->areas[0]);                        //This is done by just duplicating the 0th element.
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
            else if (testClustLeft(clust, x_size)){                               //checks 5) in the notes.
                clust->areas.push_back(clust->areas[0]);                        //This is done by just duplicating the 0th element.
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
                std::cout << clust->index << " has collided with " << other->index << std::endl;
                std::cout << "clust coord: " << clust->particles[0].pos.x << "," << clust->particles[0].pos.y << std::endl;
                std::cout << "other coord: " << other->particles[0].pos.x << "," << other->particles[0].pos.y << std::endl;

                std::cout << "clust[0] tlx = " << clust->areas[0].top_left.x << std::endl;
                std::cout << "clust[0] brx = " << clust->areas[0].bottom_right.x << std::endl;
                std::cout << "other[0] tlx = " << other->areas[0].top_left.x << std::endl;
                std::cout << "other[0] brx = " << other->areas[0].bottom_right.x << std::endl;


                std::cout << "nbr areas out = " << nbrAreasOut(clust, other, x_size, y_size) << std::endl;
                std::cout << "does the cluster we focus on cross the x-boundary? " << crossXBound(*clust, x_size) << std::endl;
                std::cout << "does the target cross the x-boundary? " << crossXBound(*other, x_size) << std::endl;
                std::cout << "does the cluster we focus on cross the y-boundary? " << crossYBound(*clust, x_size) << std::endl;
                std::cout << "does the target cross the y-boundary? " << crossYBound(*other, x_size) << std::endl;
                if (nbrAreasOut(clust, other, x_size, y_size) == 2){
                    std::cout << "is the target on the left side of the domain? " << testOtherLeft(other, x_size) << std::endl;
                    std::cout << "is the target on the right side of the domain? " << testOtherRight(other, x_size) << std::endl;
                    std::cout << "is the cluster on the left side of the domain? " << testClustLeft(clust, x_size) << std::endl;
                    std::cout << "is the cluster on the right side of the domain? " << testClustRight(clust, x_size) << std::endl;
                    std::cout << "is the target on the top side of the domain? " << testClustUp(other, y_size) << std::endl;
                    std::cout << "is the cluster on the top side of the domain? " <<testClustUp(clust, y_size) << std::endl;
                }
                std::cout << std::endl;
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
            else if (testClustUp(other, y_size)){                                 //checks 7) in the notes OTHER
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
            else if (testClustDown(other, y_size)){                               //checks 8) in the notes. OTHER
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
            if (testClustUp(clust, y_size)){                                    //checks 9) in notes. CLUST
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
            else if (testClustDown(clust, y_size)){                             //checks 10) in the notes. CLUST
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
        else if ((crossXBound(*other, x_size)) && (crossYBound(*other, y_size))
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
//                    clust->areas[1].top_left.x = other->areas[0].top_left.x;
                    clust->areas[1].top_left.x = 0;
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

velocity calc_vel(Quadtree_vel* vel_tree, Cluster* clust){
    return vel_tree->queryRange_vel(clust->particles[0].pos);
}

void takeStepTest(std::map<int, Cluster*> &clusters, Quadtree_vel &vel_tree){
    velocity tmp;
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end(); iterator++){
        if (iterator->second != nullptr){
            tmp = vel_tree.queryRange_vel(iterator->second->particles[0].pos);
            iterator->second->particles[0].pos.x += tmp.v*std::cos(tmp.theta);
            iterator->second->particles[0].pos.y += tmp.v*std::sin(tmp.theta);
        }
    }
}

void visualizeVelocity(std::vector<sf::RectangleShape> &lines,
                       std::vector<sf::CircleShape> &triangles,
                       std::vector<sf::Transform> &line_transforms,
                       std::vector<sf::Transform> &tri_transforms,
                       int vel_gen, int x_size, int y_size, double PI,
                       Quadtree_vel &tree){
    int grid_size = std::sqrt(std::pow(4, vel_gen));
    double h = std::min(x_size, y_size)/double(grid_size);
    double scale, tri_len, a;

    Point tmp_p;
    velocity tmp;
    for (int i = 0; i < grid_size; ++i){
        for (int j = 0; j < grid_size; ++j){
            tmp_p.x = j*h + h/2.0 + 0.01;
            tmp_p.y = i*h + h/2.0 + 0.01;
            tmp = tree.queryRange_vel(tmp_p);
            scale = tmp.v/2.0;
            tri_len = (h/4.0)*scale;
            a = std::sqrt(tri_len*tri_len + (tri_len/2.0)*(tri_len/2.0));

            sf::Transform transform;
            transform.rotate(180 + (180*tmp.theta)/PI, i*h + h/2.0 , j*h + h/2.0);
            sf::RectangleShape line(sf::Vector2f(h/2.0, 3));
            line.setPosition(i*h + h/2.0, j*h + h/2.0);
            lines.push_back(line);
            line_transforms.push_back(transform);

            sf::Transform tri_trans;
            tri_trans.rotate(90 + (180*tmp.theta)/PI,
                             i*h + h/2.0,
                             j*h + h/2.0).translate(i*h + h/2.0 - tri_len, j*h + h/2.0 - a/2.0);
            sf::CircleShape triangle(tri_len, 3);
            triangle.setFillColor(sf::Color::White);
            triangles.push_back(triangle);
            tri_transforms.push_back(tri_trans);

        }
    }
}

double findLifetime(eddy e){
    return e.L/std::sqrt(e.E);
}

double findStepLength(Cluster* cluster, Quadtree_vel &vel_tree, double PI,
                      double dt, double &step_dir, double diff_threshold,
                      double L_typical, double rho_air, double rho_dust,
                      double C_sphere){
    double step_len = 0;
    if (cluster != nullptr){
        double sx = 0;
        double sy = 0;
        velocity temp = vel_tree.queryRange_vel(cluster->CM);
        if (cluster->radius*L_typical < diff_threshold){
            double D = (1.38*std::pow(10.0, -23.0)*293.0*10.0)/                 //The 10 is to take cunningham correctino factor into account
                    (6.0*PI*1.75*std::pow(10.0,-5.0)*cluster->radius*L_typical);
            step_len = std::sqrt(2.0*dt*D);
            sx += step_len*std::cos(step_dir) + temp.v*std::cos(temp.theta)*dt;
            sy += step_len*std::sin(step_dir) + temp.v*std::sin(temp.theta)*dt;
            step_len = std::sqrt(std::pow(sx,2.0) + std::pow(sy, 2.0))/L_typical;
            step_dir = findDirection(sx, sy, PI);
//            std::cout << "step_len = " << step_len << " units" << std::endl;
//            std::cout << "step_len = " << step_len*L_typical << " m" << std::endl;
        }
        else if (cluster->radius*L_typical > diff_threshold){
//            std::cout << "cluster->particles.size() = " << cluster->particles.size() << std::endl;
//            std::cout << "cluster->vel.v = " << cluster->vel.v << " [m/s]" << std::endl;
//            std::cout << "cluster->vel.theta = " << cluster->vel.theta << " [rad]" << std::endl;
//            std::cout << "temp.v = " << temp.v << " [m/s]" << std::endl;
//            std::cout << "temp.theta = " << temp.theta << " [rad]" << std::endl;
//            std::cout << "v_surr_x = " << temp.v*std::cos(temp.theta) << std::endl;
//            std::cout << "v_surr_y = " << temp.v*std::sin(temp.theta) << std::endl;
//            std::cout << "v_clust_x = " << cluster->vel.v*std::cos(cluster->vel.theta) << std::endl;
//            std::cout << "v_clust_y = " << cluster->vel.v*std::sin(cluster->vel.theta) << std::endl;
//            std::cout << "cluster->radius = " << cluster->radius << std::endl;
//            std::cout << "cluster->index = " << cluster->index << std::endl;

            double m_clust = rho_dust*(4.0/3.0)*PI*std::pow(cluster->radius*L_typical,3.0);
//            std::cout << "m_clust = " << m_clust << std::endl;
            double A = PI*std::pow(cluster->radius*L_typical, 2.0);
            double Fx = 0.5*rho_air*std::pow(cluster->vel.v*
                                             std::cos(cluster->vel.theta) -
                                             temp.v*std::cos(temp.theta),2.0)*
                                             C_sphere*A;
            double Fy = 0.5*rho_air*std::pow(cluster->vel.v*
                                             std::sin(cluster->vel.theta) -
                                             temp.v*std::sin(temp.theta),2.0)*
                                             C_sphere*A;
            double a_x = Fx/m_clust;
            double a_y = Fy/m_clust;
//            std::cout << "a_x = " << a_x << " [m/s^2]" << std::endl;
//            std::cout << "a_y = " << a_y << " [m/s^2]" << std::endl;
            double D = (1.38*std::pow(10.0, -23.0)*293.0)/                      //Note no cunningham correction factor here.
                    (6.0*PI*1.75*std::pow(10.0,-5.0)*cluster->radius*L_typical);
            step_len = std::sqrt(2.0*dt*D);
            sx += step_len*std::cos(step_dir) + cluster->vel.v*std::cos(cluster->vel.theta) + (a_x*dt);
            sy += step_len*std::sin(step_dir) + cluster->vel.v*std::sin(cluster->vel.theta) + (a_y*dt);
//            sx = cluster->vel.v*std::cos(cluster->vel.theta) + (a_x*dt);
//            sy = cluster->vel.v*std::sin(cluster->vel.theta) + (a_y*dt);
            step_len = std::sqrt(std::pow(sx, 2.0) + std::pow(sy, 2.0));
            step_dir = findDirection(sx, sy, PI);
//            std::cout << "Fx = " << Fx << std::endl;
//            std::cout << "Fy = " << Fy << std::endl;
//            std::cout << "step_len = " << step_len << " [m/s]" << std::endl;
//            std::cout << "step_len = " << (step_len/L_typical)*dt << " [units/timestep]" << std::endl;
//            std::cout << "step_dir = " << step_dir << std::endl;
//            std::cout << std::endl;
            cluster->vel.v = step_len;
            cluster->vel.theta = step_dir;
            step_len = (step_len/L_typical)*dt;
        }
        else {
            std::cout << "something went wrong when calculating step length" << std::endl;
        }
    }
    return step_len;
}

double findDirection(double dx, double dy, double PI){
    double hyp = std::sqrt(dx*dx + dy*dy);
    if ((dx > 0) && (dy > 0)){
        return std::asin(dy/hyp);
    }
    else if ((dx < 0) && (dy > 0)){
        return std::acos(dx/hyp);
    }
    else if ((dx < 0) && (dy < 0)){
        if (std::abs(dx) < std::abs(dy)){
            return PI - std::asin(dy/hyp);
        }
        else {
            return 2*PI - std::acos(dx/hyp);
        }
    }
    else if ((dx > 0) && (dy < 0)){
        return std::asin(dy/hyp);
    }
}

Point FindCM(Cluster clust, double rho_carbon, double rho_dust, double r_dust,
             double r_carbon, double PI, double L_typical){
    double M = 0, m = 0;
    double sum_x = 0, sum_y = 0;
    for (auto && particle:clust.particles){
        if (particle.r_p == r_carbon){
            m = rho_carbon*(4.0/3.0)*PI*std::pow(particle.r_p*L_typical, 3.0);
            M += m;
            sum_x += m*particle.pos.x;
            sum_y += m*particle.pos.y;

        }
        else if (particle.r_p == r_dust){
            m = rho_dust*(4.0/3.0)*PI*std::pow(particle.r_p*L_typical, 3.0);
            M += m;
            sum_x += m*particle.pos.x;
            sum_y += m*particle.pos.y;
        }
    }
    return Point(sum_x/M, sum_y/M);
}
