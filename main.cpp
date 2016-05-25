#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <memory>
#include <fstream>
#include <ctime>

#include "containers.h"
#include "routines.h"

using namespace std;

int main(){
    bool visualize = false;
    bool fallout = true;
    bool varying_size = false;
    double x_size = 2200;
    double y_size = 1500;
    int nbr_particles = 20000;
    int nbr_dust = 90;
    double r_p = 1.0;
    double r_dust = 100;
    double u_0 = 0.18;
    int vel_generations = 9;

    double E_0 = std::pow(u_0, 2.0);                                            //Initial turbulent kinetic energy.
    double len = 1.001;
    double step_len = 1.0;
    double step_dir = 0;
    double PI = 3.1415926535897932384626433832;
    int col_with = -1;
    double t = 0.0;
    double p1 = 0.7;                                                            //Fraction of energy distributed to halv of daughters
    double p2 = 1.0 - p1;
    double rho_air = 1.29;                                                      //Mass density of air [kg/m^3]
    double rho_carbon = 2000.0;                                                 //Mass density of carbon particles [kg/m^3]
    double rho_dust = ((2.65+2.60+2.75+2.35+2.82+2.95+2.71+2.87+2.30)/9.0)*1000;//Mass density of dust particles [kg/m^3]
    double L_typical = 1.0*std::pow(10.0, -8.0);                                //Typical length scale of the system is set to 10 nm per unit. [m]
    double C_sphere = 0.47;
    double carbon_density = 0.0;
    double dust_density = 0.0;
    double carbon_system_size = 200000;
    double dust_system_size = 5000;
    double simulation_time = 1.0;
    double diff_threshold = 1.0*std::pow(10.0, -6.0);                           // [m]
    double dt = 1000.0/(u_0/L_typical);                                           //timestep per iteration. [s]
    double Cc = 10.0;                                                             //Cunningham correction factor
    double dyn_visc_air = 1.8*std::pow(10.0, -5.0);
    double gravity = 9.8;
    double max_lifetime = 0;
    double vt = findSettlingVelocity(rho_dust, r_dust, dyn_visc_air, gravity,
                                     1.01, L_typical, max_lifetime);
    int redist_carb = 0;
    int redist_dust = 0;
    int nbrs = nbr_particles + nbr_dust;
    int tot_fallout_carbon = 0;
    std::vector<int> fallout_sizes;
    double vel_size = 0;

    std::cout << "dt = " << dt << std::endl;
    std::cout << "vt = " << vt << std::endl;
    std::cout << "max_lifetime = " << max_lifetime << std::endl;
    if (visualize){
        r_p = 1.0;
        r_dust = 100.0;
        x_size = 2200;
        y_size = 1500;
        nbr_dust = 0;
    }
    else {
        x_size = findSystemSize(r_p, Cc, PI, dt, L_typical, vel_generations,
                                E_0, p1, simulation_time, dust_density,
                                carbon_density, nbr_particles+nbr_particles,
                                nbr_dust, nbr_particles, dust_system_size,
                                carbon_system_size, carbon_system_size,
                                dust_system_size, max_lifetime);
        y_size = x_size;
    }
    double reduced_size = x_size/L_typical;
    double red_dust_size = dust_system_size/L_typical;
    double red_carb_size = carbon_system_size/L_typical;
    std::cout << "main x_size = " << x_size << std::endl;
    std::cout << "main y_Size = " << y_size << std::endl;
    std::cout << "main carbon_size = " << carbon_system_size << std::endl;
    std::cout << "main dust_size = " << dust_system_size << std::endl;

    AABB domain = AABB(Point(0,reduced_size), Point(reduced_size,0));
    Quadtree tree(domain);

//---------------------------------Seeding RNG ---------------------------------
    int i_max = std::numeric_limits<int>::max();
    std::mt19937::result_type seed = time(0);
    auto rand_seed = std::bind(std::uniform_int_distribution<int>(0,i_max),
                       std::mt19937(seed));
    auto rand_dir = std::bind(std::uniform_real_distribution<double>(0,2*PI),
                   std::mt19937(seed));
    std::random_device rd;
    std::mt19937 g(rd());

    std::vector<Cluster> clusters;
    std::map<int, Cluster*> clusters_test;
    std::vector<Cluster> tmp_res;
    std::vector<Cluster> full_res;
    AABB search_range(Point(0,0), Point(0,0));
    distributeParticlesTest(reduced_size, reduced_size, nbr_particles, clusters,
                            rand_seed(), rand_seed(), r_p, varying_size,
                            nbr_dust, r_dust, clusters_test, rho_dust, PI,
                            rho_carbon, L_typical, red_carb_size,
                            red_dust_size, visualize, max_lifetime,
                            simulation_time);
    typedef std::map<int, Cluster*>::iterator it_type;
    vel_size = (reduced_size)/3.0;
    std::cout << "turbulence size = " << vel_size*L_typical << std::endl;
    AABB_vel vel_domain = AABB_vel(Point(0, vel_size),
                                   Point(vel_size, 0));

    eddy init_eddy(E_0, rand_dir(), vel_size);
    std::vector<double> p_list;
    p_list.push_back(p1);
    p_list.push_back(p1);
    p_list.push_back(p2);
    p_list.push_back(p2);
    Quadtree_vel vel_tree(vel_domain, init_eddy, 0, 0.0,
                              findLifetime(init_eddy, L_typical), p_list);
    vel_tree.subdivide_vel(g, rand_seed(), PI, vel_generations, L_typical);
    printTurnovertimes(vel_tree, L_typical);
    string filename = "/Users/berntsim/Documents/Master/Data/c" +
                        to_string(int(nbr_particles/1000)) + "k_d" +
                        to_string(int(nbr_dust/1000)) + "k_" +
                        to_string(int(red_dust_size/1000)) +
                        "k_" +
                         to_string(int(red_carb_size/1000))
                        +  "k.txt";
    std::ofstream out_stream;
    out_stream.open(filename);
    const clock_t begin_time = clock();
//---------------------------------Simulations----------------------------------
    if (visualize){
        return 0;
    }
    else {
        std::cout << "particles placed, start of simulations" << std::endl;
        std::cout << "Simulation time = " << simulation_time << std::endl;
        while (t < simulation_time){
//        for (int i = 0; i < 1; ++i){
            vel_tree.updateEddyTree(t, vel_generations, g, PI, rand_seed(),
                                    p_list, L_typical);
            for(it_type iterator = clusters_test.begin();
                iterator != clusters_test.end(); iterator++){
                if (iterator->second != nullptr){
                    step_dir = rand_dir();
                    len = findStepLength(iterator->second, vel_tree, PI,
                                         dt, step_dir, diff_threshold,
                                         L_typical, rho_air, C_sphere,
                                         reduced_size, vel_size);
                    for (auto&& area:iterator->second->areas){                  //The following loop prints the areas of the clusters
                        search_range = testFindSearchRange(area, len);
                        tmp_res = testBPcolCheck(search_range, tree);           //...and checks for collisions.
                        if (tmp_res.size() > 1){
                            full_res.insert(full_res.end(), tmp_res.begin(),
                                            tmp_res.end());
                        }
                        tmp_res.clear();
                    }
                    step_len = NPColCheck(iterator->second, full_res, len,
                                          step_dir, col_with, reduced_size,
                                          reduced_size, clusters_test);
                    takeSingleStep(step_dir, step_len, iterator->second,
                                   reduced_size, reduced_size);
                    if (step_len < len){
                        TestJoinClusters(iterator->second,
                                         clusters_test[col_with],
                                         clusters_test, reduced_size,
                                         reduced_size, PI, rho_carbon,
                                         rho_dust, r_dust, r_p, L_typical);
                        tree.insert(*iterator->second);
                    }
                    if (fallout){
                        fallOut(clusters_test, iterator->second, t, r_p, r_dust,
                                redist_carb, redist_dust, fallout_sizes);
                    }
                    full_res.clear();
                }
            }
            if ((fallout) &&  ((redist_carb > 0) || (redist_dust > 0))){
                redistribute(reduced_size, reduced_size, redist_carb,
                             redist_dust, rand_seed(), rand_seed(), r_p, r_dust,
                             clusters_test, PI, L_typical, simulation_time, t,
                             red_carb_size, red_dust_size, varying_size,
                             rho_carbon, rho_dust, max_lifetime, nbrs);
            }
            tree.clearQadtree();
            testPutInQuadtree(clusters_test, tree);
            writePositionsToFile(clusters_test);
            tot_fallout_carbon += redist_carb;
            out_stream << t << "    "
                       << testFindTotalCarbonOnDust(clusters_test, r_dust, r_p)
                       << "    " << largestCluster(clusters_test)
                       << "    " << meanSizeCluster(clusters_test)
                       << "    " << findTotAmountPart(clusters_test)
                       << "    " << redist_carb
                       << "    " << redist_dust
                       << "    " << tot_fallout_carbon;
            for (auto&& clust_size:fallout_sizes){
                out_stream << "    " << clust_size;
            }
            out_stream << std::endl;
            t += dt;
            redist_carb = 0;
            redist_dust = 0;
            fallout_sizes.clear();
            std::cout << "t = " << t << std::endl;
        }
        out_stream.close( );
    }
    std::cout << double( clock () - begin_time ) /  CLOCKS_PER_SEC
              << " seconds runtime" << std::endl;
    return 0;
}
