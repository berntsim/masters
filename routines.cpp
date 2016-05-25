#include "routines.h"
#include "containers.h"

double findDistance(Point a, Point b, double x_size, double y_size,
                    double L_typical){
    x_size = x_size/L_typical;
    y_size = y_size/L_typical;
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


void distributeParticlesTest(double len, double height, int nbr_particles,
                             std::vector<Cluster> &clusters,
                             std::mt19937::result_type seed_x,
                             std::mt19937::result_type seed_y, double r_p,
                             bool varying_size, int nbr_dust, double r_dust,
                             std::map<int, Cluster*> &clust_test, double rho_dust,
                             double PI, double rho_carbon, double L_typical,
                             double carbon_system_size, double dust_system_size,
                             bool visualize, double max_lifetime,
                             double simulation_time){
    double min_coord_x_c = len/6.0;
    double max_coord_x_c = (5.0*len)/6.0;
    double min_coord_y_c = height/6.0;
    double max_coord_y_c = (5.0*height)/6.0;
    double min_coord_x_d = len/6.0;
    double max_coord_x_d = (5.0*len)/6.0;
    double min_coord_y_d = height/6.0;
    double max_coord_y_d = (5.0*height)/6.0;
    if (!visualize){
        min_coord_x_c = (len/2.0 - carbon_system_size/2.0);
        max_coord_x_c = (len/2.0 + carbon_system_size/2.0);
        min_coord_y_c = (height/2.0 - carbon_system_size/2.0);
        max_coord_y_c = (height/2.0 + carbon_system_size/2.0);
        min_coord_x_d = (len/2.0 - dust_system_size/2.0);
        max_coord_x_d = (len/2.0 + dust_system_size/2.0);
        min_coord_y_d = (height/2.0 - dust_system_size/2.0);
        max_coord_y_d = (height/2.0 + dust_system_size/2.0);
    }


    std::cout << "dust x range: " << max_coord_x_d - min_coord_x_d << std::endl;
    std::cout << "dust y range: " << max_coord_y_d - min_coord_y_d << std::endl;
    std::cout << "carbon x range: " << max_coord_x_c-min_coord_x_c << std::endl;
    std::cout << "carbon y range: " << max_coord_y_c-min_coord_y_c << std::endl;

    auto coord_x_d = std::bind(std::uniform_real_distribution<double>(
                                   min_coord_x_d, max_coord_x_d),
                               std::mt19937(seed_x));
    auto coord_y_d = std::bind(std::uniform_real_distribution<double>(
                                   min_coord_y_d, max_coord_y_d),
                               std::mt19937(seed_y));


    auto coord_x = std::bind(std::uniform_real_distribution<double>(
                                 min_coord_x_c, max_coord_x_c),                 //We define a function to generate a random point in the x-dimension on the domain
                               std::mt19937(seed_x));
    auto coord_y = std::bind(std::uniform_real_distribution<double>(
                                 min_coord_y_c, max_coord_y_c),                 //similarly for the y-dimension.
                               std::mt19937(seed_y));



    auto rand_size = std::bind(std::uniform_real_distribution<double>(1.0,5.0),
                               std::mt19937(seed_x));
    auto rand_lifetime = std::bind(std::uniform_real_distribution<double>(0,
                                                                  max_lifetime),
                               std::mt19937(seed_x));
    double dist;
    int part_placed = 0;
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;
    Velocity tmp_vel;
    double m_dust = rho_dust*PI*(4.0/3.0)*std::pow(r_dust*L_typical, 3.0);
    double m_carbon = rho_carbon*PI*(4.0/3.0)*std::pow(r_p*L_typical, 3.0);

    while (int(clusters.size()) < nbr_dust){
        tmp.x = coord_x_d();
        tmp.y = coord_y_d();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clusters.size()); ++i){
            dist = findDistance(clusters[i].particles[0].pos, tmp, len, height,
                                L_typical);
            if (clusters[i].particles[0].r_p + r_dust > dist){
                occupied = true;
            }
        }
        if (!occupied){
            tmp_part.push_back(Particle(tmp, r_dust));
            aabb_tmp.bottom_right.x = tmp.x+r_dust;
            aabb_tmp.bottom_right.y = tmp.y - r_dust;
            aabb_tmp.top_left.x = tmp.x - r_dust;
            aabb_tmp.top_left.y = tmp.y + r_dust;
            areas_tmp.push_back(aabb_tmp);
            clusters.push_back(Cluster(false, areas_tmp, tmp_part, part_placed,
                                       r_dust, tmp_vel, tmp, m_dust,
                                       rand_lifetime()));
            clust_test.insert(std::make_pair(part_placed, new Cluster(false,
                                             areas_tmp, tmp_part, part_placed,
                                             r_dust, tmp_vel, tmp, m_dust,
                                             rand_lifetime())));
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;                                                       //We reset and prepare to place the next particle.
    }
    while (int(clusters.size()) < nbr_particles + nbr_dust){
        tmp.x = coord_x();
        tmp.y = coord_y();
        if (varying_size){
            r_p = rand_size();
        }
        for (int i = 0; i < int(clusters.size()); ++i){
            dist = findDistance(clusters[i].particles[0].pos, tmp, len, height,
                                L_typical);
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
            clusters.push_back(Cluster(false, areas_tmp, tmp_part, part_placed,
                                       r_p, tmp_vel, tmp, m_carbon,
                                       simulation_time));
            clust_test.insert(std::make_pair(part_placed, new Cluster(false,
                                             areas_tmp, tmp_part, part_placed,
                                             r_p, tmp_vel, tmp, m_carbon,
                                             simulation_time)));
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;                                                       //We reset and prepare to place the next particle.
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

void takeSingleStep(double step_dir, double len, Cluster* cluster, double x_size,
                    double y_size){
    if (cluster != nullptr){
        for (auto&& particle: cluster->particles){
            particle.pos.x += len*std::cos(step_dir);
            particle.pos.y += len*std::sin(step_dir);
            if (particle.pos.x < 0){
                particle.pos.x += x_size;
                std::cout << "pos = " << particle.pos.x << std::endl;
                std::cout << "We have a problem 1" << std::endl;
            }
            else if (particle.pos.x >= x_size){
                particle.pos.x -= x_size;
                std::cout << "pos = " << particle.pos.x << std::endl;
                std::cout << "We have a problem 2" << std::endl;
            }
            if (particle.pos.y < 0){
                particle.pos.y += y_size;
                std::cout << "pos = " << particle.pos.y << std::endl;
                std::cout << "We have a problem 3" << std::endl;
            }
            else if (particle.pos.y >= y_size){
                particle.pos.y -= y_size;
                std::cout << "pos = " << particle.pos.y << std::endl;
                std::cout << "We have a problem 4" << std::endl;
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

double LHit(double step_L, double step_dir, Particle one, Particle two,
            double x_size, double y_size){
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
                  double x_size, double y_size,
                  std::map<int, Cluster*> &clust_test){
    double tmp;
    double ret_val = step_len;
    if (cluster != nullptr){
        for (auto&& target : targets){
            if (clust_test.at(target.index) != nullptr){
                target = *clust_test.at(target.index);
                if (target.index != cluster->index){
                    for (auto&& particle: cluster->particles){
                        for (auto&& tar_part: target.particles){
                            tmp = LHit(step_len, step_dir, particle, tar_part, x_size,
                                       y_size);
                            if (tmp <= ret_val){
                                ret_val = tmp;
                                col_with = target.index;
                            }
                        }
                    }
                }
            }
        }
    }
    return ret_val;
}


void TestJoinClusters(Cluster* clust, Cluster* other,
                      std::map<int, Cluster*> &clusters, double x_size,
                      double y_size, double PI, double rho_carbon,
                      double rho_dust, double r_dust, double r_carbon,
                      double L_typical){
    if ((clust != nullptr) && (other != nullptr)){
        if (clust->index != other->index){
//            if (other->index >= 10500){
//                std::cout << other->index << std::endl;
//                std::cout << "a dust particle with index > 10500 has collided" << std::endl;
//            }
            updateAreas(clust, other);
            consMomentum(clust, *other, PI);
            clust->mass += other->mass;
            clust->particles.insert(clust->particles.end(),
                                    other->particles.begin(),
                                    other->particles.end());
            clust->CM = FindCM(*clust, rho_carbon, rho_dust, r_dust, r_carbon,
                               PI, L_typical);
            clust->lifetime = std::min(clust->lifetime, other->lifetime);
            double diameter = 0;
            for (int i = 0; i < int(clust->particles.size()); ++i){
                for (int j = i+1; j < int(clust->particles.size()); ++j){
                    double tmp_diameter = findDistance(clust->particles[i].pos,
                                                       clust->particles[j].pos,
                                                       x_size, y_size,
                                                       L_typical) +
                                          clust->particles[i].r_p +
                                          clust->particles[j].r_p;
                    if (tmp_diameter> diameter){
                        diameter = tmp_diameter;
                    }
                }
            }
            clust->radius = diameter/2.0;
            clusters.at(other->index) = nullptr;
        }
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

void updateAreas(Cluster* clust, Cluster* other){
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

double findLifetime(eddy e, double L_typical){
    return e.L*L_typical/std::sqrt(e.E);
}

double findStepLength(Cluster* cluster, Quadtree_vel &vel_tree, double PI,
                      double dt, double &step_dir, double diff_threshold,
                      double L_typical, double rho_air,  double C_sphere,
                      double size, double vel_size){
    double step_len = 0;
    if (cluster != nullptr){
        double sx = 0;
        double sy = 0;
        Velocity temp = vel_tree.queryRange_vel(mapToVelDomain(size, vel_size,
                                                               cluster->CM));
        if (cluster->radius*L_typical < diff_threshold){
            double D = (1.38*std::pow(10.0, -23.0)*293.0*10.0)/                 //The 10 is to take cunningham correctino factor into account
            (6.0*PI*1.75*std::pow(10.0,-5.0)*cluster->radius*L_typical);
            step_len = std::sqrt(2.0*dt*D);
            sx += step_len*std::cos(step_dir) + temp.v*std::cos(temp.theta)*dt;
            sy += step_len*std::sin(step_dir) + temp.v*std::sin(temp.theta)*dt;
            step_len = std::sqrt(std::pow(sx,2.0) + std::pow(sy, 2.0))/L_typical;
            step_dir = findDirection(sx, sy, PI);
        }
        else if (cluster->radius*L_typical >= diff_threshold){
            double A = PI*std::pow(cluster->radius*L_typical, 2.0);
            double Fx = 0.5*rho_air*std::pow(cluster->vel.v*
                                             std::cos(cluster->vel.theta) -
                                             temp.v*std::cos(temp.theta),2.0)*
                                             C_sphere*A;
            double Fy = 0.5*rho_air*std::pow(cluster->vel.v*
                                             std::sin(cluster->vel.theta) -
                                             temp.v*std::sin(temp.theta),2.0)*
                                             C_sphere*A;
            double a_x = Fx/cluster->mass;
            double a_y = Fy/cluster->mass;
            double D = (1.38*std::pow(10.0, -23.0)*293.0)/                      //Note no cunningham correction factor here.
                    (6.0*PI*1.75*std::pow(10.0,-5.0)*cluster->radius*L_typical);
            step_len = std::sqrt(2.0*dt*D);
            sx += step_len*std::cos(step_dir) +
                    cluster->vel.v*std::cos(cluster->vel.theta) + (a_x*dt);
            sy += step_len*std::sin(step_dir) +
                    cluster->vel.v*std::sin(cluster->vel.theta) + (a_y*dt);
            step_len = std::sqrt(std::pow(sx, 2.0) + std::pow(sy, 2.0));
            step_dir = findDirection(sx, sy, PI);
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
    return 0;
}

void printTurnovertimes(Quadtree_vel &vel_tree, double L_typical){
    std::cout << vel_tree.turnovertime << "  " << vel_tree.E.theta << std::endl;
    std::cout << vel_tree.ne->turnovertime << "  " << vel_tree.ne->E.theta << std::endl;
    std::cout << vel_tree.ne->ne->turnovertime << "  " << vel_tree.ne->ne->E.theta << std::endl;
    std::cout << vel_tree.ne->ne->ne->turnovertime << "  " << vel_tree.ne->ne->ne->E.theta << std::endl;
    std::cout << vel_tree.ne->ne->ne->ne->turnovertime << "  " << vel_tree.ne->ne->ne->ne->E.theta << std::endl;
    std::cout << vel_tree.ne->ne->ne->ne->ne->turnovertime << "  " << vel_tree.ne->ne->ne->ne->ne->E.theta << std::endl;
    std::cout << vel_tree.ne->ne->ne->ne->ne->ne->turnovertime << "  " << vel_tree.ne->ne->ne->ne->ne->ne->E.theta << std::endl;
    std::cout << vel_tree.ne->ne->ne->ne->ne->ne->ne->turnovertime << std::endl;
    std::cout << vel_tree.ne->ne->ne->ne->ne->ne->ne->ne->turnovertime << std::endl;

    double L_0 = (vel_tree.boundary.bottom_right.x - vel_tree.boundary.top_left.x)*L_typical;
    std::cout << "L_0 = " << L_0 << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,1.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,2.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,3.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,4.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,5.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,6.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,7.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,8.0) << " m" << std::endl;
    std::cout << L_0/std::pow(2.0,9.0) << " m" << std::endl;
}

void consMomentum(Cluster* clust, Cluster other, double PI){
    double px_c = 0, py_c = 0, px_o = 0, py_o = 0, p = 0, theta = 0;            //Momentum in x and y direction for cluster and other respectively.
    px_c = clust->mass*clust->vel.v*std::cos(clust->vel.theta);
    py_c = clust->mass*clust->vel.v*std::sin(clust->vel.theta);
    px_o = other.mass*other.vel.v*std::cos(other.vel.theta);
    py_o = other.mass*other.vel.v*std::sin(other.vel.theta);

    px_c += px_o;
    py_c += py_o;

    p = std::sqrt(px_c*px_c + py_c*py_c);
    theta = findDirection(px_c, py_c, PI);
    clust->vel.v = p/(clust->mass + other.mass);
    clust->vel.theta = theta;
}

Point FindCM(Cluster clust, double rho_carbon, double rho_dust, double r_dust,
             double r_carbon, double PI, double L_typical){
    double M = 0;
    double sum_x = 0, sum_y = 0;
    double m_carbon = rho_carbon*(4.0/3.0)*PI*std::pow(r_carbon*L_typical, 3.0);
    double m_dust = rho_dust*(4.0/3.0)*PI*std::pow(r_dust*L_typical, 3.0);
    for (auto && particle:clust.particles){
        if (particle.r_p == r_carbon){
            M += m_carbon;
            sum_x += m_carbon*particle.pos.x;
            sum_y += m_carbon*particle.pos.y;

        }
        else if (particle.r_p == r_dust){
            M += m_dust;
            sum_x += m_dust*particle.pos.x;
            sum_y += m_dust*particle.pos.y;
        }
    }
    return Point(sum_x/M, sum_y/M);
}

double findSystemSize(double r_p, double Cc, double PI, double dt,
                      double L_typical, int vel_generations, double E_0,
                      double p1, double simulation_time, double &dust_density,
                      double carbon_density, int tot_objects, int nbr_dust,
                      int nbr_carbon, double &dust_size, double &carbon_size,
                      double carbon_system_size, double dust_system_size,
                      double max_lifetime){
    double D = (1.38*std::pow(10.0, -23.0)*293.0*Cc)/
            (6.0*PI*1.75*std::pow(10.0,-5.0)*r_p*L_typical);
    double diff_step_len = std::sqrt(2.0*dt*D);
    double u_highest = 0;
    for (int i = 0; i < vel_generations; ++i){
        u_highest += std::sqrt(std::pow(p1, i)*E_0);
    }
    double max_movement = u_highest*(dt/L_typical) + diff_step_len/L_typical;   //units/timestep is the the dimensionless displacement.
    double size = 2*max_movement*(L_typical/dt)*simulation_time +
                carbon_size*L_typical;
    dust_size = 0.5;
    std::cout << "max dust = " << dust_size + 2*(u_highest*max_lifetime) << std::endl;
    std::cout << "max carbon = " << carbon_system_size*L_typical + (max_movement*L_typical/dt)*
                 simulation_time << std::endl;
    size = std::max(dust_size + 2*(u_highest*max_lifetime),
                    carbon_system_size*L_typical + (max_movement*L_typical/dt)*
                    simulation_time);
    std::cout << "u_highest = " << u_highest << std::endl;
    carbon_size = carbon_system_size*L_typical;
    std::cout << "carbon_size = " << carbon_size << std::endl;
    std::cout << "dust_size = " << dust_size << std::endl;
    std::cout << "size = " << size << std::endl;
    return size;
}

void writePositionsToFile(std::map<int, Cluster*> &clusters_test){
    std::ofstream out_stream;
    out_stream.open("/Users/berntsim/Documents/Master/Data/ParticleConfig.csv");
    out_stream.precision(10);
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters_test.begin();
        iterator != clusters_test.end(); iterator++){
        if (iterator->second != nullptr){
            for (auto&& particle:iterator->second->particles){
                out_stream << particle.pos.x << "," << particle.pos.y
                           << "," << particle.r_p << ","
                           << iterator->second->index << std::endl;
            }
        }
    }
    out_stream.close( );
}

void printNumberOfClusters(std::map<int, Cluster*> &clusters_test){
    int counter = 0;
    typedef std::map<int, Cluster*>::iterator it_type;
    for(it_type iterator = clusters_test.begin();
        iterator != clusters_test.end(); iterator++){
        if (iterator->second != nullptr){
            counter++;
        }
    }
    std::cout << "counter = " << counter << std::endl;
}

int testFindTotalCarbonOnDust(std::map<int, Cluster*> &clusters_test,
                              double r_dust, double r_carbon){
    typedef std::map<int, Cluster*>::iterator it_type;
    int nbr_carbon = 0;
    bool isDust = false;
    for(it_type iterator = clusters_test.begin();
        iterator != clusters_test.end(); iterator++){
        if (iterator->second != nullptr){
            for (auto&& particle:iterator->second->particles){
                if (particle.r_p == r_dust){
                    isDust = true;
                }
            }
            for (auto&& particle:iterator->second->particles){
                if ((particle.r_p == r_carbon) && (isDust)){
                    nbr_carbon++;
                }
            }
            isDust = false;
        }
    }
    return nbr_carbon;
}

int largestCluster(std::map<int, Cluster*> &clusters_test){
    typedef std::map<int, Cluster*>::iterator it_type;
    int largest_cluster = 0;
    int current_size = 0;
    for(it_type iterator = clusters_test.begin();
        iterator != clusters_test.end(); iterator++){
        if (iterator->second != nullptr){
            current_size = iterator->second->particles.size();
            if (current_size > largest_cluster){
                largest_cluster = current_size;
            }
        }
    }
    return largest_cluster;
}

double meanSizeCluster(std::map<int, Cluster*> &clusters_test){
    typedef std::map<int, Cluster*>::iterator it_type;
    int tot = 0;
    int nbr_clusters = 0;
    for(it_type iterator = clusters_test.begin();
        iterator != clusters_test.end(); iterator++){
        if (iterator->second != nullptr){
            nbr_clusters++;
            tot += iterator->second->particles.size();
        }
    }
    return double(tot)/double(nbr_clusters);
}


double findSettlingVelocity(double rho_dust, double r_dust, double dyn_visc_air,
                            double g, double Cc, double L_typical,
                            double &max_lifetime){
    double vt = (1.0/18.0)*(std::pow(2.0*r_dust*L_typical, 2.0)*rho_dust*g*Cc)/
            dyn_visc_air;
    max_lifetime = (r_dust*L_typical)/(vt);
    return vt;
}

void fallOut(std::map<int, Cluster*> &clusters_test, Cluster* clust, double t,
             double r_carbon, double r_dust, int &redist_carbon,
             int &redist_dust, std::vector<int> &fallout_sizes){
    if (clust != nullptr){
        if (t > clust->lifetime){
//            std::cout << "cluster number " << clust->index
//                      << " had a lifetime of " << clust->lifetime
//                      << " seconds and has now fallen out."
//                      << std::endl;
            redist_carbon += countParticleOccurence(clust->particles,
                                                    r_carbon);
            redist_dust += countParticleOccurence(clust->particles, r_dust);
            fallout_sizes.push_back(redist_carbon);
            clusters_test.at(clust->index) = nullptr;
        }
    }
}

int countParticleOccurence(std::vector<Particle> particles, double r_particle){
    int counter = 0;
    for (auto&& particle:particles){
        if (particle.r_p == r_particle){
            ++counter;
        }
    }
    return counter;
}


void redistribute(double x_size, double y_size, int nbr_carbon, int nbr_dust,
                  std::mt19937::result_type seed_x,
                  std::mt19937::result_type seed_y, double r_carbon,
                  double r_dust, std::map<int, Cluster*> &clust_test,
                  double PI, double L_typical, double simulation_time, double t,
                  double size_carbon, double size_dust, bool varying_size,
                  double rho_carbon, double rho_dust, double max_lifetime,
                  int &org_size){
    double min_coord_x_c = (x_size/2.0 - size_carbon/2.0);
    double max_coord_x_c = (x_size/2.0 + size_carbon/2.0);
    double min_coord_y_c = (y_size/2.0 - size_carbon/2.0);
    double max_coord_y_c = (y_size/2.0 + size_carbon/2.0);
    double min_coord_x_d = (x_size/2.0 - size_dust/2.0);
    double max_coord_x_d = (x_size/2.0 + size_dust/2.0);
    double min_coord_y_d = (y_size/2.0 - size_dust/2.0);
    double max_coord_y_d = (y_size/2.0 + size_dust/2.0);

    auto coord_x_d = std::bind(std::uniform_real_distribution<double>(
                                   min_coord_x_d, max_coord_x_d),
                               std::mt19937(seed_x));
    auto coord_y_d = std::bind(std::uniform_real_distribution<double>(
                                   min_coord_y_d, max_coord_y_d),
                               std::mt19937(seed_y));
    auto coord_x_c = std::bind(std::uniform_real_distribution<double>(
                                 min_coord_x_c, max_coord_x_c),
                               std::mt19937(seed_x));
    auto coord_y_c = std::bind(std::uniform_real_distribution<double>(
                                 min_coord_y_c, max_coord_y_c),
                               std::mt19937(seed_y));
    auto rand_size = std::bind(std::uniform_real_distribution<double>(1.0,5.0),
                               std::mt19937(seed_x));
    double dist = 0;
    int part_placed = 0;
    Point tmp;
    bool occupied = false;
    std::vector<Particle> tmp_part;
    AABB aabb_tmp = AABB(Point(), Point());
    std::vector<AABB> areas_tmp;
    Velocity tmp_vel;
    int col_with = -1;
    double m_dust = rho_dust*(4.0/3.0)*PI*std::pow(r_dust*L_typical, 3.0);
    double m_carbon = rho_carbon*(4.0/3.0)*PI*std::pow(r_carbon*L_typical, 3.0);
//    std::cout << "nbr_carbon = " << nbr_carbon << std::endl;
//    std::cout << "nbr_dust = " << nbr_dust << std::endl;
//    std::cout << "reached start of distribution" << std::endl;
    while (part_placed < nbr_dust){
        tmp.x = coord_x_d();
        tmp.y = coord_y_d();
        for (int i = 0; i < org_size; ++i){
//            std::cout << "clust test size dust: " << clust_test.size() << std::endl;
//            std::cout << i << std::endl;
            if (clust_test.at(i) != nullptr){
//                std::cout << "the clust_test.at(i) != nullptr test passed dust" << std::endl;
//                std::cout << clust_test.at(i)->particles.size() << std::endl;
                for (auto&& particle:clust_test.at(i)->particles){
//                    std::cout << "we made it into the particles loop dust" << std::endl;
                    dist = findDistance(particle.pos, tmp, x_size, y_size,
                                        L_typical);
//                    std::cout << "we checked distance dust" << std::endl;
                    if (particle.r_p + r_dust > dist){
                        occupied = true;
                        col_with = i;
                    }
                }
//                std::cout << "we exited the particle loop dust" << std::endl;
            }
        }
//        std::cout << "finished dust for loop with one .at()" << std::endl;
        if (!occupied){
            tmp_part.push_back(Particle(tmp, r_dust));
            aabb_tmp.bottom_right.x = tmp.x + r_dust;
            aabb_tmp.bottom_right.y = tmp.y - r_dust;
            aabb_tmp.top_left.x = tmp.x - r_dust;
            aabb_tmp.top_left.y = tmp.y + r_dust;
            areas_tmp.push_back(aabb_tmp);
            clust_test.insert(std::make_pair(part_placed + org_size,
                                             new Cluster(false,
                                             areas_tmp, tmp_part, part_placed + org_size,
                                             r_dust, tmp_vel, tmp, m_dust,
                                             t + max_lifetime)));
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
//            std::cout << "finished !occupied dust" << std::endl;
        }
        else {
            tmp_part.push_back(Particle(tmp, r_dust));
            aabb_tmp.bottom_right.x = tmp.x + r_dust;
            aabb_tmp.bottom_right.y = tmp.y - r_dust;
            aabb_tmp.top_left.x = tmp.x - r_dust;
            aabb_tmp.top_left.y = tmp.y + r_dust;
            areas_tmp.push_back(aabb_tmp);
            clust_test.insert(std::make_pair(part_placed + org_size,
                                             new Cluster(false,
                                             areas_tmp, tmp_part, part_placed + org_size,
                                             r_dust, tmp_vel, tmp, m_dust,
                                             t + max_lifetime)));
            TestJoinClusters(clust_test.at(part_placed + org_size),
                             clust_test.at(col_with), clust_test,
                             x_size/L_typical, y_size/L_typical, PI, rho_carbon,
                             rho_dust, r_dust, r_carbon, L_typical);
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
        }
//        std::cout << "finished else dust" << std::endl;
        occupied = false;
    }
    while (part_placed < nbr_carbon + nbr_dust){
        tmp.x = coord_x_c();
        tmp.y = coord_y_c();
        for (int i = 0; i < org_size; ++i){
//            std::cout << "clust test size carbon: " << clust_test.size() << std::endl;
//            std::cout << i << std::endl;
            if (clust_test.at(i) != nullptr){
//                std::cout << "the clust_test.at(i) != nullptr test passed carbon" << std::endl;
//                std::cout << clust_test.at(i)->particles.size() << std::endl;
                for (auto&& particle:clust_test.at(i)->particles){
//                    std::cout << "we made it into the particles loop carbon" << std::endl;
                    dist = findDistance(particle.pos, tmp, x_size, y_size,
                                        L_typical);
//                    std::cout << "we checked distance carbon" << std::endl;
                    if (particle.r_p + r_carbon > dist){
                        occupied = true;
                        col_with = i;
                    }
                }
//                std::cout << "we exited the particle loop carbon" << std::endl;
            }
        }
//        std::cout << "finished carbon for loop" << std::endl;
        if (!occupied){
            tmp_part.push_back(Particle(tmp, r_carbon));
            aabb_tmp.bottom_right.x = tmp.x + r_carbon;
            aabb_tmp.bottom_right.y = tmp.y - r_carbon;
            aabb_tmp.top_left.x = tmp.x - r_carbon;
            aabb_tmp.top_left.y = tmp.y + r_carbon;
            areas_tmp.push_back(aabb_tmp);
            clust_test.insert(std::make_pair(part_placed + org_size,
                                             new Cluster(false,
                                             areas_tmp, tmp_part, part_placed + org_size,
                                             r_carbon, tmp_vel, tmp, m_carbon,
                                             t + simulation_time)));
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
//            std::cout << "finished !occupied carbon" << std::endl;
        }
        else {
            tmp_part.push_back(Particle(tmp, r_carbon));
            aabb_tmp.bottom_right.x = tmp.x + r_carbon;
            aabb_tmp.bottom_right.y = tmp.y - r_carbon;
            aabb_tmp.top_left.x = tmp.x - r_carbon;
            aabb_tmp.top_left.y = tmp.y + r_carbon;
            areas_tmp.push_back(aabb_tmp);
            clust_test.insert(std::make_pair(part_placed + org_size,
                                             new Cluster(false,
                                             areas_tmp, tmp_part, part_placed + org_size,
                                             r_carbon, tmp_vel, tmp, m_carbon,
                                             t + simulation_time)));
            TestJoinClusters(clust_test.at(part_placed + org_size),
                             clust_test.at(col_with), clust_test, x_size,
                             y_size, PI, rho_carbon, rho_dust, r_dust, r_carbon,
                             L_typical);
            part_placed++;
            tmp_part.clear();
            areas_tmp.clear();
        }
        occupied = false;
//        std::cout << "finished occupied carbon" << std::endl;
    }
    org_size += nbr_carbon;
    org_size += nbr_dust;
//    std::cout << "org_size = " << org_size << std::endl;
}

int findTotAmountPart(std::map<int, Cluster*> &clusters_test){
    typedef std::map<int, Cluster*>::iterator it_type;
    int tot = 0;
    for(it_type iterator = clusters_test.begin();
        iterator != clusters_test.end(); iterator++){
        if (iterator->second != nullptr){
            tot += iterator->second->particles.size();
        }
    }
    return tot;
}


Point mapToVelDomain(double size, double vel_size, Point CM){
//    std::cout << "CM.x = " << CM.x << std::endl;
//    std::cout << "CM.y = " << CM.y << std::endl;
//    std::cout << "size = " << size << std::endl;
    double vel_min = (size - vel_size)/2.0;
    double vel_max = size - vel_min;
//    std::cout << "vel_min = " << vel_min << std::endl;
//    std::cout << "vel_max = " << vel_max << std::endl;
    Point ret = CM;
    ret.x = std::fmod(CM.x, vel_size);
    ret.y = std::fmod(CM.y, vel_size);


//    if (CM.x < vel_min){
//        ret.x = CM.x + vel_size;
//    }
//    else if (CM.x >= vel_max){
//        ret.x = CM.x - vel_size;
//    }
//    if (CM.y < vel_min){
//        ret.y = CM.y + vel_size;
//    }
//    else if (CM.y >= vel_max){
//        ret.y = CM.y - vel_size;
//    }
    if ((ret.x == 0) || (ret.y == 0)){
        std::cout << "OH SHIT" << std::endl;
    }
    return ret;
}
