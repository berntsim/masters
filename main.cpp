#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <random>
#include <map>
#include <memory>

#include "containers.h"
#include "routines.h"

using namespace std;

int main(){
    bool visualize = true;
    bool varying_size = false;
    double x_size = 2200;
    double y_size = 1500;
    int nbr_particles = 10000;
    double gamma = 1.0;
    double D_0 = 10;
    double r_p = 1.0;
    double len = 1.001;
    double step_len = 1.0;
    int iteration_length = 5;
    double step_dir = 0;
    double PI = 3.1415926535897932384626433832;
    int col_with = -1;
    double dt = 0.1;

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

    std::vector<Cluster> clusters;
    std::map<int, Cluster*> clusters_test;
    std::vector<Cluster> tmp_res;
    std::vector<Cluster> full_res;
    std::vector<AABB> search_ranges;
    std::vector<sf::VertexArray> lines_vec;
    AABB search_range(Point(0,0), Point(0,0));
    distributeParticlesTest(x_size, y_size, nbr_particles, clusters, rand_seed(),
                        rand_seed(), r_p, varying_size, clusters_test);
    typedef std::map<int, Cluster*>::iterator it_type;

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
            for(it_type iterator = clusters_test.begin();
                iterator != clusters_test.end(); iterator++){
                if (iterator->second != nullptr){
                    step_dir = rand_dir();
//                    len = calcStepLen(iterator->second, gamma, PI, D_0, dt);
                    for (auto&& area:iterator->second->areas){                  //The following loop prints the areas of the clusters
                        search_range = testFindSearchRange(area, step_len);
//                        sf::VertexArray lines = printSearchRange(search_range);
//                        lines_vec.push_back(lines);
                        tmp_res = testBPcolCheck(search_range, tree);           //...and checks for collisions.
                        if (tmp_res.size() > 1){
                            full_res.insert(full_res.end(), tmp_res.begin(),
                                            tmp_res.end());
                        }
                        tmp_res.clear();
                    }
                    step_len = NPColCheck(iterator->second, full_res, len,
                                          step_dir, col_with, x_size, y_size);
                    takeSingleStep(step_dir, step_len, iterator->second,
                                   x_size, y_size);
                    if (step_len < len){
                        TestJoinClusters(iterator->second,
                                         clusters_test[col_with],
                                         clusters_test, x_size, y_size);
                        step_dir = -PI/2.0;
                    }
//                    if (full_res.size() > 1){                                       //If collisions occur, the clusters are drawn red
//                        addToDraw(collided, full_res, col_fill_color);
//                    }
                    full_res.clear();
                    if (iterator->second != nullptr){
                        splitAreas(*iterator->second, x_size, y_size);
                    }
                }
            }

        tree.clearQadtree();
        testPutInQuadtree(clusters_test, tree);
        testAddToDraw(to_draw, clusters_test, fill_color);
        for (auto&& draw:to_draw){
            window.draw(draw);
        }
//        for (auto&& lines:lines_vec){
//            window.draw(lines);
//        }
//        for (auto&& draw:collided){
//            window.draw(draw);
//        }

        window.display();
        to_draw.clear();
        lines_vec.clear();
        collided.clear();
        }
    }
    else {
//        for (int i = 0; i < iteration_length; ++i){
//            std::cout << i << std::endl;

//            for(it_type iterator = clusters_test.begin();
//                iterator != clusters_test.end(); iterator++){
//                if (iterator->second != nullptr){
//                    step_dir = rand_dir();
//                    for (auto&& area:iterator->second->areas){                  //The following loop prints the areas of the clusters
//                        search_range = testFindSearchRange(area, 3*step_len);
//                        tmp_res = testBPcolCheck(search_range, tree);           //...and checks for collisions.
//                        if (tmp_res.size() > 1){
//                            full_res.insert(full_res.end(), tmp_res.begin(),
//                                            tmp_res.end());
//                        }
//                        tmp_res.clear();
//                    }
//                    step_len = NPColCheck(iterator->second, full_res, len,
//                                          step_dir, col_with, x_size, y_size);
//                    takeSingleStep(step_dir, step_len, iterator->second,
//                                   x_size, y_size);
//                    if (step_len < len){
//                        TestJoinClusters(iterator->second,
//                                         clusters_test[col_with],
//                                         clusters_test, x_size, y_size);
//                        step_dir = -PI/2.0;
//                    }
//                    full_res.clear();
//                    if (iterator->second != nullptr){
//                        splitAreas(*iterator->second, x_size, y_size);
//                    }
//                }
//            }
//            tree.clearQadtree();
//            testPutInQuadtree(clusters_test, tree);
//        }
        double test_len = calcStepLen(clusters_test.at(0), gamma, PI, D_0, dt);
        std::cout << "r_p = " << clusters_test.at(0)->particles[0].r_p << std::endl;
        std::cout << "A = " << PI*pow(clusters_test.at(0)->particles[0].r_p, 2) << std::endl;
        std::cout << "step_len = " << test_len << std::endl;
    }
    return 0;
}
