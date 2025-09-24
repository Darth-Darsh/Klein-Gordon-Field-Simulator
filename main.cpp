#include<iostream>
#include<vector>
#include<map>
#include<complex>
#include<string>
#include<iomanip>
#include<Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct> 
#include <unsupported/Eigen/MatrixFunctions> 
using namespace std::complex_literals;
using namespace std;

using complex_d = complex<double>;

class KleinGordonSimulator{
    public:
    KleinGordonSimulator(int n, double J, double h) : n(n), J(J), h(h), dim(1LL << n) {
        GATES_["I"] = Eigen::Matrix2cd::Identity();
        GATES_["X"] = (Eigen::Matrix2cd() << 0, 1, 1, 0).finished();
        GATES_["Y"] = (Eigen::Matrix2cd() << 0, -1i, 1i, 0).finished();
        GATES_["Z"] = (Eigen::Matrix2cd() << 1, 0, 0, -1).finished();

        state_vector = Eigen::VectorXcd::Zero(dim);
        state_vector(0) = 1.0;

        build_hamiltonian();
    }

    void initial_particle(int site){
        if(site >=n || site < 0){
            throw std::runtime_error("Error: site index out of range.");
        }
        state_vector = Eigen::VectorXcd::Zero(dim);
        long long index = 1LL << site;
        state_vector(index) = 1.0;
    }

    void evolve(double time){
        Eigen::MatrixXcd U = (-1i * H * time).exp();
        state_vector = U * state_vector;
    }

    vector<double> measure(){
        vector<double> probabilities(n, 0.0);
        for(long long i = 0; i < n; ++i){
            Eigen::MatrixXcd n_op = get_operator_for_site("Z", i);
            complex_d expectation_z = (state_vector.adjoint() * n_op * state_vector)(0);
            
            probabilities[i] = 0.5*(1.0 - expectation_z.real());
        }
        return probabilities;
    }
private:
    int n;
    double J;
    double h;
    long long dim;

    Eigen::VectorXcd state_vector;
    Eigen::MatrixXcd H;
    std::map<string, Eigen::Matrix2cd> GATES_;

    void build_hamiltonian() {
        H = Eigen::MatrixXcd::Zero(dim, dim);

        for (int j = 0; j < n; ++j) {
            H += h * get_operator_for_site("Z", j);
        }

        for (int j = 0; j < n - 1; ++j) {
            Eigen::MatrixXcd XX = get_operator_for_two_sites("X", j, "X", j + 1);
            Eigen::MatrixXcd YY = get_operator_for_two_sites("Y", j, "Y", j + 1);
            H += -J * (XX + YY);
        }
    }

    Eigen::MatrixXcd get_operator_for_site(const std::string& gate, int site) {
        std::vector<Eigen::Matrix2cd> op_list(n, GATES_["I"]);
        op_list[site] = GATES_[gate];
        
        Eigen::MatrixXcd full_matrix = op_list[0];
        for (int i = 1; i < n; ++i) {
            full_matrix = Eigen::kroneckerProduct(full_matrix, op_list[i]).eval();
        }
        return full_matrix;
    }

    Eigen::MatrixXcd get_operator_for_two_sites(const std::string& g1, int s1, const std::string& g2, int s2) {
        std::vector<Eigen::Matrix2cd> op_list(n, GATES_["I"]);
        op_list[s1] = GATES_[g1];
        op_list[s2] = GATES_[g2];

        Eigen::MatrixXcd full_matrix = op_list[0];
        for (int i = 1; i < n; ++i) {
            full_matrix = Eigen::kroneckerProduct(full_matrix, op_list[i]).eval();
        }
        return full_matrix;
    }
};

void stream_data_for_plotting(double time, const std::vector<double>& probs) {
    std::cout << time; 
    for(const auto& p : probs) {
        std::cout << "," << p; 
    }
    std::cout << std::endl; 
}


int main() {
    const int N_SITES = 7;
    const double J_COUPLING = 1.0;   
    const double H_FIELD = 0.5;
    const double TOTAL_TIME = 100;
    const double TIME_STEP = 1;

    // std::cout << "--- Simulating Klein-Gordon Field on a " << N_SITES << "-site lattice ---" << std::endl;
    // std::cout << "Hopping J = " << J_COUPLING << ", Mass h = " << H_FIELD << std::endl;
    // std::cout << "---------------------------------------------------------" << std::endl;

    KleinGordonSimulator sim(N_SITES, J_COUPLING, H_FIELD);

    int middle_site = N_SITES / 2;
    sim.initial_particle(middle_site);
    
    double current_time = 0.0;
    while (current_time <= TOTAL_TIME) {
        auto probs = sim.measure();
        stream_data_for_plotting(current_time, probs);
        
        sim.evolve(TIME_STEP);
        current_time += TIME_STEP;
    }

    return 0;
}

