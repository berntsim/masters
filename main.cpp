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
    double x_size = 2500;
    double y_size = 1500;
    int nbr_particles = 100000;
    double r_p = 1.0;
    double len = 1.0;
    double step_len = 1.0;
    double step_dir = 0;
    double PI = 3.1415926535897932384626433832;
    int col_with = -1;

    sf::Color fill_color = sf::Color::Green;
    sf::Color fill_color_collision = sf::Color::Red;
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


    AABB search_range = AABB(Point(x_size/2.0 - 100, y_size/2.0 + 100),
                             Point(x_size/2.0 + 100, y_size/2.0 - 100));
    AABB search_range_clust = search_range;


    std::vector<Cluster> clusters;
    std::map<int, Cluster*> clusters_test;
    std::vector<Cluster> tmp_res;
    distributeParticlesTest(x_size, y_size, nbr_particles, clusters, rand_seed(),
                        rand_seed(), r_p, varying_size, clusters_test);
//    std::cout << "clusters_test.size() = " << clusters_test.size() << std::endl;
//    std::cout << "clusters.size() = " << clusters.size() << std::endl;


    typedef std::map<int, Cluster*>::iterator it_type;
//    for(it_type iterator = clusters_test.begin(); iterator != clusters_test.end(); iterator++) {
//        std::cout << iterator->second->area.top_left.y << std::endl;
//    }



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
            for(it_type iterator = clusters_test.begin(); iterator != clusters_test.end(); iterator++){
                if (iterator->second != nullptr){
                    step_dir = rand_dir();
                    findSearchRange(search_range_clust, iterator->second, len);
                    tmp_res = BPcolCheck(iterator->second, search_range_clust, tree);
                    step_len = NPColCheck(iterator->second, tmp_res, len, step_dir, col_with);
                    takeSingleStep(step_dir, step_len, iterator->second, x_size, y_size);
                    if (step_len < len){
                        TestJoinClusters(iterator->second, clusters_test[col_with], clusters_test, col_with);
                    }
//                addToDrawRec(collided_rec, tmp_res, fill_color_collision);
                }
                tmp_res.clear();
            }
            tree.clearQadtree();
            testPutInQuadtree(clusters_test, tree);
            testAddToDraw(to_draw, clusters_test, fill_color);
//            testAddToDrawRec(to_draw_rec, clusters_test, fill_color);

            for (auto&& draw:to_draw){
                window.draw(draw);
            }

//            for (auto&& draw: collided_rec){
//                window.draw(draw);
//            }

            window.display();
            collided.clear();
            to_draw.clear();
            to_draw_rec.clear();
            collided_rec.clear();
//            std::cout << "clusters_test.size() = " << clusters_test.size() << std::endl;
        }
    }
    else {
//        for (int i = 0; i < 2; ++i){
//            for (auto&& cluster: clusters){
//                step_dir = rand_dir();
//                findSearchRange(search_range_clust, cluster, len);
//                tmp_res = BPcolCheck(cluster, search_range_clust, tree);
////                step_len = NPColCheck(cluster, tmp_res, x_size, y_size, len, step_dir);
//                tmp_res.clear();
//                takeSingleStep(step_dir, step_len, cluster, x_size, y_size);
//            }

//            tree.clearQadtree();
//            putInQuadtree(clusters, tree);
//            std::cout << "Iteration number: " << i << std::endl;
//        }
//    }
    }

    for(it_type iterator = clusters_test.begin(); iterator != clusters_test.end(); iterator++){
        delete iterator->second;
    }
    return 0;
}
