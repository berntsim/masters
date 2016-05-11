#ifndef ROUTINES_H
#define ROUTINES_H

#include <iostream>
#include <vector>
#include <SFML/Graphics.hpp>
#include <random>
#include <map>
#include <memory>

#include "containers.h"
#include "vel_field.h"


double findDistance(Point a, Point b, double x_size, double y_size);

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

//void addToDrawRec(std::vector<sf::RectangleShape> &to_draw,
//               std::vector<Cluster> clusters,
//               sf::Color fill_color);

sf::VertexArray printSearchRange(AABB range);

void takeSingleStep(double step_dir, double len, Cluster* cluster, int x_size,
                    int y_size);

std::vector<Cluster> BPcolCheck(Cluster* cluster,
                                std::vector<AABB> search_ranges,
                                Quadtree &tree);

double LHit(double step_L, double step_dir, Particle one, Particle two,
            int x_size, int y_size);

double NPColCheckOrg(Cluster* cluster, std::vector<Cluster> targets,
                  double step_len, double step_dir, int &col_with,
                  int x_size, int y_size);

double NPColCheck(Cluster* cluster, std::vector<Cluster> targets,
                  double &step_len, double &step_dir, int &col_with,
                  int x_size, int y_size, Quadtree_vel &vel_tree, double dt,
                  double L_typical, double rho_air, double rho_dust,
                  double C_sphere, double PI, bool &will_collide);













void testPutInQuadtree(std::map<int, Cluster*> clusters, Quadtree &qtree);

void testAddToDrawRec(std::vector<sf::RectangleShape> &to_draw,
                      std::map<int, Cluster*> clusters,
                      sf::Color fill_color);

void testAddToDraw(std::vector<sf::CircleShape> &to_draw,
               std::map<int, Cluster*> clusters,
               sf::Color fill_color);

void TestJoinClusters(Cluster* clust, Cluster* other,
                      std::map<int, Cluster*> &clusters,
                      int x_size, int y_size);

void printInformation(Cluster* clust, int col_with);

bool checkPointInAABB(Point pt, AABB box);

bool crossXBound(Cluster clust, int x_size);

bool crossYBound(Cluster clust, int y_size);

void splitAreas(Cluster &clust, int x_size, int y_size);

AABB testFindSearchRange(AABB area, double step_len);

std::vector<Cluster> testBPcolCheck(AABB search_range, Quadtree &tree);

int nbrAreasOut(Cluster* one, Cluster* other, int x_size, int y_size);

bool otherLeftTwo(Cluster* clust, Cluster* other);

bool clustLeftTwo(Cluster* clust, Cluster* other);

bool otherRightTwo(Cluster* clust, Cluster* other);

bool clustRightTwo(Cluster* clust, Cluster* other);

bool otherUpTwo(Cluster* clust, Cluster* other);

bool otherDownTwo(Cluster* clust, Cluster* other);

bool otherLeftFour(Cluster* clust, Cluster* other);

bool otherRightFour(Cluster* clust, Cluster* other);

bool otherDownFour(Cluster* clust, Cluster* other);

bool otherUpFour(Cluster* clust, Cluster* other);

void updateAreas(Cluster* clust, Cluster* other, int x_size, int y_size);

velocity calc_vel(Quadtree_vel *vel_tree, Cluster *clust);

void takeStepTest(std::map<int, Cluster*> &clusters, Quadtree_vel &vel_tree);

void visualizeVelocity(std::vector<sf::RectangleShape> &lines,
                       std::vector<sf::CircleShape> &triangles,
                       std::vector<sf::Transform> &line_transforms,
                       std::vector<sf::Transform> &tri_transforms,
                       int vel_gen, int x_size, int y_size, double PI,
                       Quadtree_vel &tree);

double findLifetime(eddy e);

void distributeDust(int len, int height, int nbr_particles,
                    std::mt19937::result_type seed_x, bool varying_size,
                    std::mt19937::result_type seed_y, double r_p,
                    std::map<int, Cluster*> &clusters);

double findStepLength(Cluster* cluster, Quadtree_vel &vel_tree, double PI,
                      double dt, double &step_dir, double diff_threshold,
                      double L_typical, double rho_air, double rho_dust,
                      double C_sphere);

double findDirection(double dx, double dy, double PI);

#endif // ROUTINES_H
