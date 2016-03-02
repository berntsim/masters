#ifndef ROUTINES_H
#define ROUTINES_H

#include <iostream>
#include <vector>
#include <SFML/Graphics.hpp>
#include <random>
#include <map>
#include <memory>

#include "containers.h"


//void distributeParticles(int len, int height, int nbr_particles,
//                            std::vector<Cluster> &clusters,
//                            std::mt19937::result_type seed_x,
//                            std::mt19937::result_type seed_y, double r_p,
//                            bool varying_size);

void distributeParticlesTest(int len, int height, int nbr_particles,
                            std::vector<Cluster> &clusters,
                            std::mt19937::result_type seed_x,
                            std::mt19937::result_type seed_y, double r_p,
                            bool varying_size, std::map<int, Cluster*> &clust_test);

void addToDraw(std::vector<sf::CircleShape> &to_draw,
               std::vector<Cluster> clusters,
               sf::Color fill_color);

void addToDrawRec(std::vector<sf::RectangleShape> &to_draw,
               std::vector<Cluster> clusters,
               sf::Color fill_color);

sf::VertexArray printSearchRange(AABB range);

void putInQuadtree(std::vector<Cluster> clusters, Quadtree &qtree);

void takeStep(std::mt19937::result_type seed, double len,
              std::vector<Cluster> &clusters,int x_size, int y_size);

void takeSingleStep(double step_dir, double len, Cluster* cluster, int x_size,
                    int y_size);

void findSearchRange(AABB &range, Cluster* cluster, double step_len);

std::vector<Cluster> BPcolCheck(Cluster* cluster, AABB search_range,
                                Quadtree &tree);

double LHit(double step_L, double step_dir, Particle one, Particle two);

double NPColCheck(Cluster* cluster, std::vector<Cluster> targets,
                  double step_len, double step_dir, int &col_with);

void joinClusters(Cluster &clust, Cluster &other, std::vector<Cluster> &clusters);













void testPutInQuadtree(std::map<int, Cluster*> clusters, Quadtree &qtree);

void testAddToDrawRec(std::vector<sf::RectangleShape> &to_draw,
                      std::map<int, Cluster*> clusters,
                      sf::Color fill_color);

void testAddToDraw(std::vector<sf::CircleShape> &to_draw,
               std::map<int, Cluster*> clusters,
               sf::Color fill_color);

void TestJoinClusters(Cluster* clust, Cluster* other, std::map<int, Cluster*> &clusters, int col_with);

#endif // ROUTINES_H
