#ifndef ROUTINES_H
#define ROUTINES_H

#include <iostream>
#include <vector>
#include <random>
#include <map>
#include <memory>
#include <fstream>

#include "containers.h"


//void distributeParticles(int len, int height, int nbr_particles,
//                            std::vector<Cluster> &clusters,
//                            std::mt19937::result_type seed_x,
//                            std::mt19937::result_type seed_y, double r_p,
//                            bool varying_size);

void distributeParticlesTest(double len, double height, int nbr_particles,
                             std::vector<Cluster> &clusters,
                             std::mt19937::result_type seed_x,
                             std::mt19937::result_type seed_y, double r_p,
                             bool varying_size, int nbr_dust, double r_dust,
                             std::map<int, Cluster*> &clust_test, double rho_dust,
                             double PI, double rho_carbon, double L_typical,
                             double carbon_system_size, double dust_system_size,
                             bool visualize, double max_lifetime,
                             double simulation_time);

void testPutInQuadtree(std::map<int, Cluster*> clusters, Quadtree &qtree);

void takeSingleStep(double step_dir, double len, Cluster* cluster, double x_size,
                    double y_size);

double LHit(double step_L, double step_dir, Particle one, Particle two,
            double x_size, double y_size);

double NPColCheck(Cluster* cluster, std::vector<Cluster> targets,
                  double step_len, double step_dir, int &col_with,
                  double x_size, double y_size,
                  std::map<int, Cluster*> &clust_test);

void TestJoinClusters(Cluster* clust, Cluster* other,
                      std::map<int, Cluster*> &clusters, double x_size,
                      double y_size, double PI, double rho_carbon,
                      double rho_dust, double r_dust, double r_carbon,
                      double L_typical);

bool crossXBound(Cluster clust, int x_size);

bool crossYBound(Cluster clust, int y_size);

AABB testFindSearchRange(AABB area, double step_len);

std::vector<Cluster> testBPcolCheck(AABB search_range, Quadtree &tree);

int nbrAreasOut(Cluster* one, Cluster* other, int x_size, int y_size);

void updateAreas(Cluster* clust, Cluster* other, double x_size, double y_size,
                 int cai, int oai);

double findLifetime(eddy e, double L_typical);

double findStepLength(Cluster* cluster, Quadtree_vel &vel_tree, double PI,
                      double dt, double &step_dir, double diff_threshold,
                      double L_typical, double rho_air, double C_sphere,
                      double size, double vel_size);

double findDirection(double dx, double dy, double PI);

void printTurnovertimes(Quadtree_vel &vel_tree, double L_typical);

void consMomentum(Cluster* clust, Cluster other, double PI);

Point FindCM(Cluster clust, double rho_carbon, double rho_dust, double r_dust,
             double r_carbon, double PI, double L_typical);

double findSystemSize(double r_p, double Cc, double PI, double dt,
                      double L_typical, int vel_generations, double E_0,
                      double p1, double simulation_time, double &dust_density,
                      double carbon_density, int tot_objects, int nbr_dust,
                      int nbr_carbon, double &dust_size, double &carbon_size,
                      double carbon_system_size, double dust_system_size,
                      double max_lifetime);

void writePositionsToFile(std::map<int, Cluster*> &clusters_test);

void printNumberOfClusters(std::map<int, Cluster*> &clusters_test);

int testFindTotalCarbonOnDust(std::map<int, Cluster*> &clusters_test,
                              double r_dust, double r_carbon);

int largestCluster(std::map<int, Cluster*> &clusters_test);

double meanSizeCluster(std::map<int, Cluster*> &clusters_test);

double findSettlingVelocity(double rho_dust, double r_dust, double dyn_visc_air,
                            double g, double Cc, double L_typical,
                            double &max_lifetime);

void fallOut(std::map<int, Cluster*> &clusters_test, Cluster* clust, double t,
             double r_carbon, double r_dust, int &redist_carbon,
             int &redist_dust, std::vector<int> &fallout_sizes);

int countParticleOccurence(std::vector<Particle> particles, double r_particle);

void redistribute(double x_size, double y_size, int nbr_carbon, int nbr_dust,
                  std::mt19937::result_type seed_x,
                  std::mt19937::result_type seed_y, double r_carbon,
                  double r_dust, std::map<int, Cluster*> &clust_test,
                  double PI, double L_typical, double simulation_time, double t,
                  double size_carbon, double size_dust, bool varying_size,
                  double rho_carbon, double rho_dust, double max_lifetime,
                  int &org_size);

int findTotAmountPart(std::map<int, Cluster*> &clusters_test);

Point mapToVelDomain(double size, double vel_size, Point CM);

void insertResult(std::vector<Cluster> tmp_res, std::vector<Cluster> &full_res,
                  int idx);

int determineRegion(Cluster clust, double x_size, double y_size);

void mapTo(int clust_idx, int other_idx, Cluster* other, double x_size,
           double y_size);

int intersectingIdx(std::vector<AABB> areas_clust,
                    std::vector<AABB> areas_other, int &clust_index);

#endif // ROUTINES_H
