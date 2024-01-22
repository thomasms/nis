#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>
#include <ENDFtk.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/endf/core.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/cram.hpp"

using namespace nis;

int main(int argc, char **argv){

    // cli args
    std::string filename = "default";
    int nuclide_zai = 0;
    double dt = 0.;

    CLI::App app{"NIS"};
    app.add_option("-f,--file", filename, "An ENDF file")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-n,--nuclide", nuclide_zai, "Nuclide to examine")
      ->required();
    app.add_option("-t,--dt", dt, "Time to examine")
      ->required();

    CLI11_PARSE(app, argc, argv);

    // construct a decay matrix - straight from ENDF file(s)
    auto tape = njoy::ENDFtk::tree::fromFile(filename);
    const Eigen::MatrixX<double> dense_decay_mat = endf::construct_decay_matrix_from_tape(tape);

    // set initial vector (inventory)
    const int vec_size = static_cast<int>(dense_decay_mat.cols());
    Eigen::VectorX<double> n0 = Eigen::VectorX<double>::Zero(vec_size);
    const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);
    if(zai_map.count(nuclide_zai) == 1){
        const auto zai_index = zai_map.at(nuclide_zai);
        n0(zai_index) = 1e20;
    }

    // we should actually use a sparse matrix here really
    // but for testing we make a dense mat first and convert to sparse
    const Eigen::MatrixX<double> dense_decay_mat_ld = dense_decay_mat.cast<double>();
    const Eigen::SparseMatrix<double> A = dense_decay_mat_ld.sparseView();
    const Eigen::SparseMatrix<core::cram48::Coeff_t> A_complex = A.transpose()*dt*core::cram48::Coeff_t(1,0);
    const Eigen::SparseMatrix<double> I = Eigen::MatrixX<double>::Identity(dense_decay_mat.cols(), dense_decay_mat.rows()).sparseView();

    // solver
    Eigen::SparseLU<Eigen::SparseMatrix<core::cram48::Coeff_t>, Eigen::COLAMDOrdering<int> > lu_solver;
    Eigen::SparseMatrix<core::cram48::Coeff_t> A_prime;
    Eigen::VectorX<core::cram48::Coeff_t> x(dense_decay_mat.cols());
    core::cram48::Coeff_t alpha, theta;

    Eigen::VectorX<double> nt = n0*core::cram48::ALPHA_0;
    Eigen::VectorX<core::cram48::Coeff_t> nt_c = nt*core::cram48::Coeff_t(1,0);

    // See page 4 of https://www.osti.gov/servlets/purl/1362192
    // for implementation
    for(int i=0;i<core::cram48::N_COEFF;++i){
      alpha = core::cram48::ALPHA[i];
      theta = core::cram48::THETA[i];
      A_prime = A_complex - theta*I;
      A_prime.makeCompressed();
      lu_solver.compute(A_prime);
      if(lu_solver.info() != Eigen::ComputationInfo::Success) {
        // decomposition failed
        std::cout << "Solver decomposition failed!\n";
        return -1;
      }
      x = lu_solver.solve(nt_c);
      nt += 2.0*(alpha*x).real();
      nt_c = nt*core::cram48::Coeff_t(1,0);
    }

    // output
    for(int i=0;i<vec_size;++i){
      double val = nt[i];
      int zai = (i == vec_size-1) ? -1 : data::ALL_ZAIS[i];
      if(val != 0){
        std::cout << nuclide_zai << " -> " << zai << " " << i << " (" << val << ")\n";
      }
    }

    return 0;
}
