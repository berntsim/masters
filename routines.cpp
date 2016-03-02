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

//void distributeParticles(int len, int height, int nbr_particles,
//                            std::vector<Cluster> &clusters,
//                            std::mt19937::result_type seed_x,
//                            std::mt19937::result_type seed_y, double r_p,
//                            bool varying_size){
//    auto coord_x = std::bind(std::uniform_real_distribution<double>(0,len),     //We define a function to generate a random point in the x-dimension on the domain
//                               std::mt19937(seed_x));
//    auto coord_y = std::bind(std::uniform_real_distribution<double>(0,height),  //similarly for the y-dimension.
//                               std::mt19937(seed_y));
//    auto rand_size = std::bind(std::uniform_real_distribution<float>(1.0,5.0),
//                               std::mt19937(seed_x));
//    double dist;
//    int part_placed = 0;
//    Point tmp;
//    bool occupied = false;
//    std::vector<Particle> tmp_part;
//    AABB aabb_tmp = AABB(Point(), Point());


//    while (int(clusters.size()) < nbr_particles){
//        tmp.x = coord_x();
//        tmp.y = coord_y();
//        if (varying_size){
//            r_p = rand_size();
//        }
//        for (int i = 0; i < int(clusters.size()); ++i){
//            dist = findDistance(clusters[i].particles[0].pos, tmp, len, height);
//            if (clusters[i].particles[0].r_p + r_p > dist){
//                occupied = true;
//            }
//        }
//        if (!occupied){
//            tmp_part.push_back(Particle(tmp, r_p));
//            aabb_tmp.bottom_right.x = tmp.x+r_p;
//            aabb_tmp.bottom_right.y = tmp.y - r_p;
//            aabb_tmp.top_left.x = tmp.x - r_p;
//            aabb_tmp.top_left.y = tmp.y + r_p;
//            clusters.push_back(Cluster(aabb_tmp, tmp_part, part_placed));
//            part_placed++;
//            tmp_part.clear();
//        }
//        occupied = false;                                                       //We reset and prepare to place the next particle.
//    }
//}

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
//            Cluster tmp_clust = Cluster(aabb_tmp, tmp_part, part_placed);
            clusters.push_back(Cluster(false, aabb_tmp, tmp_part, part_placed));
            clust_test.insert(std::make_pair(part_placed, new Cluster(false, aabb_tmp, tmp_part, part_placed)));
            part_placed++;
            tmp_part.clear();
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

void addToDrawRec(std::vector<sf::RectangleShape> &to_draw,
               std::vector<Cluster> clusters,
               sf::Color fill_color){
    for (auto&& clust: clusters){
        sf::RectangleShape shape(sf::Vector2f((clust.area.bottom_right.x - clust.area.top_left.x),
                                              (clust.area.top_left.y - clust.area.bottom_right.y)));
        shape.setPosition(clust.area.top_left.x, clust.area.top_left.y-(clust.area.top_left.y - clust.area.bottom_right.y));
        shape.setFillColor(fill_color);
        to_draw.push_back(shape);
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
        cluster.area.top_left.x += len*std::cos(dir);
        cluster.area.bottom_right.x += len*std::cos(dir);
        cluster.area.top_left.y += len*std::sin(dir);
        cluster.area.bottom_right.y += len*std::sin(dir);
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
        cluster->area.top_left.x += len*std::cos(step_dir);
        cluster->area.bottom_right.x += len*std::cos(step_dir);
        cluster->area.top_left.y += len*std::sin(step_dir);
        cluster->area.bottom_right.y += len*std::sin(step_dir);
    }
}

void findSearchRange(AABB &range, Cluster* cluster, double step_len){
    if (cluster != nullptr){
        range.bottom_right.x = cluster->area.bottom_right.x + step_len;
        range.bottom_right.y = cluster->area.bottom_right.y - step_len;
        range.top_left.x = cluster->area.top_left.x - step_len;
        range.top_left.y = cluster->area.top_left.y + step_len;
    }
    //Should change this to include periodic boundary conditions, i. e. should allow for up to 4 possible search ranges,
    //in the case of 4 possible hashes.
}

std::vector<Cluster> BPcolCheck(Cluster* cluster, AABB search_range,
                                Quadtree &tree){
    std::vector<Cluster> ret_vec;
    std::vector<Cluster> tmp_res;
    tmp_res = tree.queryRange(search_range);
    if (cluster != nullptr){
        if (tmp_res.size() > 1){
            for (auto&& res: tmp_res){
                if (!(res.area.top_left.x == cluster->area.top_left.x) &&
                    !(res.area.top_left.y == cluster->area.top_left.y)){
                    ret_vec.push_back(res);
                }
            }
        }
    }
    return ret_vec;
}

double LHit(double step_L, double step_dir, Particle one, Particle two){
    double d_p = one.r_p + two.r_p;
    double a = 1.0;
    double b = 2*(std::cos(step_dir)*(one.pos.x - two.pos.x) +
                  std::sin(step_dir)*(one.pos.y - two.pos.y));
    double c = (two.pos.x - one.pos.x)*(two.pos.x - one.pos.x) +
               (two.pos.y - one.pos.y)*(two.pos.y - one.pos.y) - d_p*d_p;
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


double NPColCheck(Cluster* cluster, std::vector<Cluster> targets,
                  double step_len, double step_dir, int &col_with){
    double tmp;
    double ret_val = step_len;
    if (cluster != nullptr){
        for (auto&& target : targets){
            for (auto&& particle: cluster->particles){
                for (auto&& tar_part: target.particles){
                    tmp = LHit(step_len, step_dir, particle, tar_part);
                    if (tmp < ret_val){
                        ret_val = tmp;
                        col_with = target.index;
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
    if (other.area.top_left.x < clust.area.top_left.x){
        clust.area.top_left.y = other.area.top_left.y;
    }
    if (other.area.top_left.y > clust.area.top_left.y){
        clust.area.top_left.y = other.area.top_left.y;
    }
    if (other.area.bottom_right.x > clust.area.bottom_right.x){
        clust.area.bottom_right.x = other.area.bottom_right.x;
    }
    if (other.area.bottom_right.y < clust.area.bottom_right.y){
        clust.area.bottom_right.y = other.area.bottom_right.y;
    }
//    clusters.erase(clusters.begin() + other.index);
}













void testAddToDrawRec(std::vector<sf::RectangleShape> &to_draw,
                      std::map<int, Cluster*> clusters,
                      sf::Color fill_color){
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end(); iterator++){
        sf::RectangleShape shape(sf::Vector2f((iterator->second->area.bottom_right.x - iterator->second->area.top_left.x),
                                              (iterator->second->area.top_left.y - iterator->second->area.bottom_right.y)));
        shape.setPosition(iterator->second->area.top_left.x, iterator->second->area.top_left.y-(iterator->second->area.top_left.y - iterator->second->area.bottom_right.y));
        shape.setFillColor(fill_color);
        to_draw.push_back(shape);
    }
}

void testAddToDraw(std::vector<sf::CircleShape> &to_draw,
               std::map<int, Cluster*> clusters,
               sf::Color fill_color){
    sf::CircleShape shape(1.0);
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters.begin(); iterator != clusters.end(); iterator++){
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

void TestJoinClusters(Cluster* clust, Cluster* other, std::map<int, Cluster*> &clusters,
                      int col_with){
//    std::cout << "col_with = " << col_with << std::endl;
    if ((clust != nullptr) && (other != nullptr)){
        if (clust->index != other->index){
            clust->particles.insert(clust->particles.end(),
                                    other->particles.begin(),
                                    other->particles.end());
            if (other->area.top_left.x < clust->area.top_left.x){
                clust->area.top_left.x = other->area.top_left.x;
            }
            if (other->area.top_left.y > clust->area.top_left.y){
                clust->area.top_left.y = other->area.top_left.y;
            }
            if (other->area.bottom_right.x > clust->area.bottom_right.x){
                clust->area.bottom_right.x = other->area.bottom_right.x;
            }
            if (other->area.bottom_right.y < clust->area.bottom_right.y){
                clust->area.bottom_right.y = other->area.bottom_right.y;
            }
            clusters.at(other->index) = nullptr;
//    clusters.erase(other.index);
//    clusters.at(other.index)->joined = false;
        }
    }
}
