#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <random>
#include <map>
#include <memory>
#include "time.h"

#include "containers.h"
#include "routines.h"
#include "vel_field.h"

using namespace std;

int main(){
    clock_t t1,t2;
    t1=clock();
//--------------------Governing parameters for simulation-----------------------
    bool test_environment = true;                                               //The case for testing stuff.
    bool visualize = false;
    bool varying_size = false;
    bool var_dust_size = false;
    bool velocity_field = true;
    int iteration_length = 7;                                                   //If there are no visualization going on, then this parameter controls runtime

    double L_typical = 1.0*std::pow(10.0,-8.0);                                 //The typical length of the system, used to non-dimensionlize the problem.
//    double E_0 = std::pow(0.18, 2.0);                                               //The initial energy of the largest eddy.
//    double L_0 = 0.01/L_typical;                                                //The physical length of the system.
    double E_0 = std::pow(0.5, 2.0);
    double L_0 = 0.01;
    double system_size = 1400;                                                  //The systme size with respect to the displayed window.
    int nbr_particles = 1000;                                                    //Number of BC particles, who will later join together as clusters.
    int nbr_dust_particles = 3;                                                //Number of mineral dust particles.
    double r_p = 2.0*std::pow(10.0,-8.0)/L_typical;                             //Radius of BC monomers. By definition unity in these simulations (L_typical)
    double r_dust = (1.0*std::pow(10.0,-6.0))/L_typical;                        //Radius of dust particles.
    double D_0 = 10;                                                            //Diffusivity. TO FIX!!!
    double rho_air = 1.29;                                                      //Mass density of air.
    double rho_dust = 2.0;                                                      //Mass density of dust. TO FIX!!!
    double C_sphere = 0.47;                                                     //drag coefficient for spherical shape.
    int vel_generations = 8;                                                    //Essentialy sets the "resolution" of the turbulent velocity field.
    double diff_threshold = 2.0*std::pow(10.0, -7.0);                           //Sets limit for when one takes diffusion into account.

//--------------------------Constants for the simulations-----------------------
    double x_size = system_size;                                                //Defines the size of the system for visualization.
    double y_size = system_size;
    double gamma = 1.0;                                                         //Not used atm. Defines diffusion-mass relation.
    double len = 1.00;                                                          //Default step len for diffusion. Recalculated at later stage.
    double step_len = 1.0;
    double step_dir = 0;                                                        //just initializing the step direction, drawn/calculated later.
    double PI = 3.1415926535897932384626433832;                                 //Pi.
    int col_with = -1;                                                          //Initializing the ID og collision target to be invalid at the start.
    double t = 0;                                                               //Starting time. Physical time during simulation.
    double dt = (system_size*L_typical)/(std::sqrt(E_0));                       //Time step.
    double p1 = 0.7;                                                            //Defines the model used for turbulence. The multiplicative increment
    double p2 = 1.0-p1;                                                         //is taken from the original article.

//----------------------------Containers initialization-------------------------
    sf::Color fill_color = sf::Color::Green;
    sf::Color col_fill_color = sf::Color::Red;
    std::vector<sf::CircleShape> to_draw;
    std::vector<sf::RectangleShape> to_draw_rec;
    std::vector<sf::CircleShape> collided;
    std::vector<sf::RectangleShape> collided_rec;
    AABB domain = AABB(Point(0,y_size), Point(x_size,0));
    Quadtree tree(domain);


//---------------------------------Seeding RNG ---------------------------------
    int i_max = std::numeric_limits<int>::max();
    std::mt19937::result_type seed = time(0);
    auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
                       std::mt19937(seed));
    auto rand_dir = std::bind(std::uniform_real_distribution<float>(0,2*PI),
                   std::mt19937(seed));
    std::random_device rd;
    std::mt19937 g(rd());

//-----------------------------Declaring containers ----------------------------
    std::vector<Cluster> clusters;
    std::map<int, Cluster*> clusters_test;
    std::vector<Cluster> tmp_res;
    std::vector<Cluster> full_res;
    std::vector<AABB> search_ranges;
    std::vector<sf::VertexArray> lines_vec;
    AABB search_range(Point(0,0), Point(0,0));
    distributeDust(x_size, y_size, nbr_dust_particles, rand_seed(),
                   var_dust_size, rand_seed(), r_dust, clusters_test);
    distributeParticlesTest(x_size, y_size, nbr_particles, clusters, rand_seed(),
                        rand_seed(), r_p, varying_size, clusters_test);
    typedef std::map<int, Cluster*>::iterator it_type;
    AABB_vel vel_domain = AABB_vel(Point(0,y_size), Point(x_size,0));
    eddy init_eddy(E_0, rand_dir(), L_0);
    std::vector<float> p_list;
    p_list.push_back(p1);
    p_list.push_back(p1);
    p_list.push_back(p2);
    p_list.push_back(p2);
    Quadtree_vel vel_tree(vel_domain, init_eddy, 0, 0.0,
                          findLifetime(init_eddy), p_list);
    vel_tree.subdivide_vel(g, rand_seed(), PI, vel_generations);
    std::vector<sf::RectangleShape> lines;
    std::vector<sf::Transform> line_transforms;
    std::vector<sf::CircleShape> triangles;
    std::vector<sf::Transform> tri_transforms;
//-----------------------------Visualisation ----------------------------
    if (test_environment){                                                      //This is used only for testing purposes.
        clusters_test.clear();
        double bcx = x_size - 500, bcy = y_size - 100;
        double bc2x = x_size - 400, bc2y = y_size - 100;
        double dustx = x_size - 50, dusty = y_size - 50;
        velocity vel;
        Particle bc_part(Point(bcx, bcy), r_p);
        Particle bc2_part(Point(bc2x, bc2y), r_p);
        Particle dust_part(Point(dustx, dusty), r_dust);
        std::vector<Particle> bc_vec;
        std::vector<Particle> bc2_vec;
        std::vector<Particle> dust_vec;
        bc_vec.push_back(bc_part);
        bc2_vec.push_back(bc2_part);
        dust_vec.push_back(dust_part);
        std::vector<AABB> bc_areas;
        std::vector<AABB> bc2_areas;
        std::vector<AABB> dust_areas;
        AABB bc_area(Point(bcx-r_p, bcy+r_p), Point(bcx+r_p, bcy-r_p));
        AABB bc2_area(Point(bc2x-r_p, bc2y+r_p), Point(bc2x+r_p, bc2y-r_p));
        AABB dust_area(Point(dustx-r_dust, dusty+r_dust), Point(dustx+r_dust, dusty-r_dust));
        bc_areas.push_back(bc_area);
        bc2_areas.push_back(bc2_area);
        dust_areas.push_back(dust_area);
//        Cluster bc_clust(false, bc_areas, bc_particles, 0, r_p, vel);
//        Cluster dust_clust(false, dust_areas, dust_particles, 1, r_dust, vel);
        clusters_test.insert(std::make_pair(0, new Cluster(false, bc_areas, bc_vec, 0, r_p*L_typical, vel)));
        clusters_test.insert(std::make_pair(1, new Cluster(false, dust_areas, dust_vec, 1, r_dust*L_typical, vel)));
        clusters_test.insert(std::make_pair(2, new Cluster(false, bc2_areas, bc2_vec, 3, r_p*L_typical, vel)));

        sf::RenderWindow window(sf::VideoMode(x_size, y_size), "DLCA");
        window.setVerticalSyncEnabled(true);
        while (window.isOpen()){
            sf::Event event;
            while (window.pollEvent(event)){
                if (event.type == sf::Event::Closed)
                    window.close();
            }
            window.clear(sf::Color::Black);
            tree.clearQadtree();
            testPutInQuadtree(clusters_test, tree);
            for(it_type iterator = clusters_test.begin();
                iterator != clusters_test.end(); iterator++){
                if (iterator->second != nullptr){
//                    step_dir = rand_dir();
                    if (iterator->second->radius < diff_threshold){
                        step_dir = 0;
                        len = 1.0;
                    }
                    else if (iterator->second->radius > diff_threshold){
                        step_dir = PI;
                        len = 0;
                    }
//                    len = findStepLength(iterator->second, vel_tree, PI,
//                                              dt, step_dir, diff_threshold,
//                                              L_typical, rho_air, rho_dust,
//                                              C_sphere);
                    for (auto&& area:iterator->second->areas){                  //The following loop prints the areas of the clusters
                        search_range = testFindSearchRange(area, len);
                        sf::VertexArray lines = printSearchRange(search_range);
                        lines_vec.push_back(lines);
                        tmp_res = testBPcolCheck(search_range, tree);           //...and checks for collisions.
                        if (tmp_res.size() > 1){
                            full_res.insert(full_res.end(), tmp_res.begin(),
                                            tmp_res.end());
                        }
                        tmp_res.clear();
                    }
                    std::cout << clusters_test.at(2)->areas.size() << std::endl;
                    step_len = NPColCheck(iterator->second, full_res, len,
                                          step_dir, col_with, x_size, y_size,
                                          vel_tree, dt, L_typical, rho_air,
                                          rho_dust, C_sphere, PI);
                    takeSingleStep(step_dir, step_len, iterator->second, x_size,
                                   y_size);
                    if (step_len < len){
                        std::cout << "hi" << std::endl;
                        std::cout << "col with = " << col_with << std::endl;
                        TestJoinClusters(iterator->second,
                                         clusters_test[col_with],
                                         clusters_test, x_size, y_size);
                    }
                    full_res.clear();
                    splitAreas(*iterator->second, x_size, y_size);
                }
            }
            testAddToDraw(to_draw, clusters_test, fill_color);
            for (auto&& draw:to_draw){
                window.draw(draw);
            }
            for (auto&& lines:lines_vec){
                window.draw(lines);
            }
            to_draw.clear();
            lines_vec.clear();
            window.display();
        }
        return 0;
    }

    if (visualize){
        sf::RenderWindow window(sf::VideoMode(x_size, y_size), "DLCA");
        window.setVerticalSyncEnabled(true);
        while (window.isOpen()){
            sf::Event event;
            while (window.pollEvent(event)){
                if (event.type == sf::Event::Closed)
                    window.close();
            }
            window.clear(sf::Color::Black);
            tree.clearQadtree();
            testPutInQuadtree(clusters_test, tree);
            for(it_type iterator = clusters_test.begin();
                iterator != clusters_test.end(); iterator++){
                if (iterator->second != nullptr){
                    step_dir = rand_dir();
                    len = findStepLength(iterator->second, vel_tree, PI,
                                              dt, step_dir, diff_threshold,
                                              L_typical, rho_air, rho_dust,
                                              C_sphere);
//                    len = 100.0;
                    for (auto&& area:iterator->second->areas){                  //The following loop prints the areas of the clusters
                        search_range = testFindSearchRange(area, len);
                        sf::VertexArray lines = printSearchRange(search_range);
                        lines_vec.push_back(lines);
                        tmp_res = testBPcolCheck(search_range, tree);           //...and checks for collisions.
                        if (tmp_res.size() > 1){
                            full_res.insert(full_res.end(), tmp_res.begin(),
                                            tmp_res.end());
                        }
                        tmp_res.clear();
                    }
                    step_len = NPColCheck(iterator->second, full_res, len,
                                          step_dir, col_with, x_size, y_size,
                                          vel_tree, dt, L_typical, rho_air,
                                          rho_dust, C_sphere, PI);
                    takeSingleStep(step_dir, step_len, iterator->second, x_size,
                                   y_size);
                    if (step_len < len){
                        TestJoinClusters(iterator->second,
                                         clusters_test[col_with],
                                         clusters_test, x_size, y_size);
                    }
                    full_res.clear();
                    splitAreas(*iterator->second, x_size, y_size);
                }
            }
            vel_tree.updateEddyTree(t, vel_generations, g, PI, rand_seed(),
                                    p_list);
//            visualizeVelocity(lines, triangles, line_transforms, tri_transforms,
//                              vel_generations, x_size, y_size, PI, vel_tree);
//            for (int i = 0; i < int(lines.size()); ++i){
//                window.draw(lines[i], line_transforms[i]);
//                window.draw(triangles[i], tri_transforms[i]);
//            }
            testAddToDraw(to_draw, clusters_test, fill_color);
            for (auto&& draw:to_draw){
                window.draw(draw);
            }
            for (auto&& lines:lines_vec){
                window.draw(lines);
            }
            window.display();
            to_draw.clear();
            triangles.clear();
            tri_transforms.clear();
            lines.clear();
            line_transforms.clear();
            lines_vec.clear();
            t += dt;
//            std::cout << "t = " << t << std::endl;
        }
        return 0;
    }
    else if (velocity_field){
//        std::cout << "level 0: " << vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y << std::endl;
//        std::cout << "level 1: " << vel_tree.ne->boundary.top_left.y-vel_tree.ne->boundary.bottom_right.y << std::endl;
//        std::cout << "level 2: " << vel_tree.ne->ne->boundary.top_left.y-vel_tree.ne->ne->boundary.bottom_right.y << std::endl;
//        std::cout << "level 3: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,3.0) << std::endl;
//        std::cout << "level 4: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,4.0) << std::endl;
//        std::cout << "level 5: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,5.0) << std::endl;
//        std::cout << "level 6: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,6.0) << std::endl;
//        std::cout << "level 7: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,7.0) << std::endl;
//        std::cout << "level 8: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,8.0) << std::endl;
//        std::cout << "level 9: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,9.0) << std::endl;
//        std::cout << "level 10: " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::pow(2.0,10.0) << std::endl;
//        std::cout << "t0 = " << (vel_tree.boundary.top_left.y-vel_tree.boundary.bottom_right.y)/std::sqrt(vel_tree.E.E) << std::endl;
//        std::cout << "t1 = " << (vel_tree.ne->boundary.top_left.y-vel_tree.ne->boundary.bottom_right.y)/std::sqrt(vel_tree.ne->E.E) << std::endl;
//        std::cout << "t2 = " << (vel_tree.ne->ne->boundary.top_left.y-vel_tree.ne->ne->boundary.bottom_right.y)/std::sqrt(vel_tree.ne->ne->E.E) << std::endl;
//        std::cout << "t3 = " << (vel_tree.ne->ne->ne->boundary.top_left.y-vel_tree.ne->ne->ne->boundary.bottom_right.y)/std::sqrt(vel_tree.ne->ne->ne->E.E) << std::endl;
//        std::cout << "t4 = " << (vel_tree.ne->ne->ne->ne->boundary.top_left.y-vel_tree.ne->ne->ne->ne->boundary.bottom_right.y)/std::sqrt(vel_tree.ne->ne->ne->ne->E.E) << std::endl;
//        std::cout << "t5 = " << (vel_tree.ne->ne->ne->ne->ne->boundary.top_left.y-vel_tree.ne->ne->ne->ne->ne->boundary.bottom_right.y)/std::sqrt(vel_tree.ne->ne->ne->ne->ne->E.E) << std::endl;
//        std::cout << "findLifetime_0: " << findLifetime(vel_tree.E) << std::endl;
//        std::cout << vel_tree.E.L << std::endl;
//        std::cout << vel_tree.E.E << std::endl;
//        std::cout << "findLifetime_1: " << findLifetime(vel_tree.ne->E) << std::endl;
//        std::cout << "findLifetime_2: " << findLifetime(vel_tree.ne->ne->E) << std::endl;
//        std::cout << "findLifetime_3: " << findLifetime(vel_tree.ne->ne->ne->E) << std::endl;
//        std::cout << "findLifetime_4: " << findLifetime(vel_tree.ne->ne->ne->ne->E) << std::endl;
//        std::cout << "findLifetime_5: " << findLifetime(vel_tree.ne->ne->ne->ne->ne->E) << std::endl;
        std::cout << "direction at lvl 0: " << vel_tree.E.theta << std::endl;
        std::cout << "direction at lvl 1: " << vel_tree.ne->E.theta << std::endl;
        std::cout << "direction at lvl 2: " << vel_tree.ne->ne->E.theta << std::endl;
        std::cout << "direction at lvl 3: " << vel_tree.ne->ne->ne->E.theta << std::endl;
        std::cout << vel_tree.queryRange_vel(Point(1399, 1399)).theta << std::endl;
        calc_vel(&vel_tree, clusters_test.at(0));
//        std::cout << calc_vel(&vel_tree, clusters_test.at(0)).v << std::endl;
        sf::ContextSettings settings;
        settings.antialiasingLevel = 8;

        sf::RenderWindow window(sf::VideoMode(x_size, y_size), "Velocity field",
                                sf::Style::Default, settings);
        window.setVerticalSyncEnabled(true);
        while (window.isOpen()){
            sf::Event event;
            while (window.pollEvent(event)){
                if (event.type == sf::Event::Closed)
                    window.close();
            }
            window.clear(sf::Color::Black);
//            initial_eddy.theta = rand_dir();
//            vel_tree.clearQuadtree_vel();
            vel_tree.updateEddyTree(t, vel_generations, g, PI, rand_seed(),
                                    p_list);
            takeStepTest(clusters_test, vel_tree);

            visualizeVelocity(lines, triangles, line_transforms, tri_transforms,
                              vel_generations, x_size, y_size, PI, vel_tree);
            for (int i = 0; i < int(lines.size()); ++i){
                window.draw(lines[i], line_transforms[i]);
                window.draw(triangles[i], tri_transforms[i]);
            }

            testAddToDraw(to_draw, clusters_test, fill_color);
            for (auto&& draw:to_draw){
                window.draw(draw);
            }

            window.display();
            to_draw.clear();
            triangles.clear();
            tri_transforms.clear();
            lines.clear();
            line_transforms.clear();
            t += dt;
//            std::cout << "t = " << t << std::endl;
        }
    }
    else {
        for (int i = 0; i < iteration_length; ++i){
            std::cout << i << std::endl;
            for(it_type iterator = clusters_test.begin();
                iterator != clusters_test.end(); iterator++){
                if (iterator->second != nullptr){
                    step_dir = rand_dir();
                    for (auto&& area:iterator->second->areas){                  //The following loop prints the areas of the clusters
                        search_range = testFindSearchRange(area, 3*step_len);
                        tmp_res = testBPcolCheck(search_range, tree);           //...and checks for collisions.
                        if (tmp_res.size() > 1){
                            full_res.insert(full_res.end(), tmp_res.begin(),
                                            tmp_res.end());
                        }
                        tmp_res.clear();
                    }
                    step_len = NPColCheck(iterator->second, full_res, len,
                                          step_dir, col_with, x_size, y_size,
                                          vel_tree, dt, L_typical, rho_air,
                                          rho_dust, C_sphere, PI);
                    takeSingleStep(step_dir, step_len, iterator->second,
                                   x_size, y_size);
                    if (step_len < len){
                        TestJoinClusters(iterator->second,
                                         clusters_test[col_with],
                                         clusters_test, x_size, y_size);
                        step_dir = -PI/2.0;
                    }
                    full_res.clear();
                    if (iterator->second != nullptr){
                        splitAreas(*iterator->second, x_size, y_size);
                    }
                }
            }
            tree.clearQadtree();
            testPutInQuadtree(clusters_test, tree);
        }

        double test_len = calcStepLen(clusters_test.at(0), gamma, PI, D_0, dt);
        std::cout << "r_p = " << clusters_test.at(0)->particles[0].r_p << std::endl;
        std::cout << "A = " << PI*pow(clusters_test.at(0)->particles[0].r_p, 2) << std::endl;
        std::cout << "step_len = " << test_len << std::endl;
    }
    t2=clock();
    std::cout << "This run took " << (t2-t1)/CLOCKS_PER_SEC << " seconds" << std::endl;
    return 0;
}
