#ifndef ENDF_REACTION_HPP
#define ENDF_REACTION_HPP

#include <functional>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include <ENDFtk.hpp>

#include <Eigen/Core>

#include "nis/spectrum.hpp"
#include "nis/zai.hpp"

namespace nis {
    namespace endf {

        // (primary, {secondaries}) - same layout as DecayProductT
        // a primary of -1 means the parent goes to the sink
        using ReactionProductT = std::tuple<int, std::vector<int>>;

        constexpr int PROTON = 10010;
        constexpr int DEUTERON = 10020;
        constexpr int TRITON = 10030;
        constexpr int HELION = 20030;

        // 1 barn = 1e-24 cm^2, the flux is assumed to be in n/cm^2/s
        constexpr double BARN_TO_CM2 = 1e-24;

        // MT=18 (total fission) has no yields yet, so we send it to the sink
        constexpr int FISSION_MT = 18;

        // Change in Z and A from the target to the residual nucleus,
        // for an incident neutron, plus the light charged particles emitted.
        // Emitted neutrons and gammas are not tracked.
        struct ReactionChange {
            int dz;
            int da;
            std::vector<int> secondaries;
        };

        // We only include the partial reactions that change the nuclide.
        // Sums (MT=1,3,4,5,27,101), elastic (MT=2), inelastic levels (MT=51-91)
        // and fission chances (MT=19-21,38) are left out on purpose,
        // otherwise we would double count or the nuclide does not change.
        // NOTE: all residuals are in the ground state, isomer production needs MF=9/10
        const std::map<int, ReactionChange> REACTION_MAP = {
            {16,  {0,  -1, {}}},                                  // (n,2n)
            {17,  {0,  -2, {}}},                                  // (n,3n)
            {22,  {-2, -4, {zai::ALPHA}}},                        // (n,na)
            {23,  {-6, -12, {zai::ALPHA, zai::ALPHA, zai::ALPHA}}}, // (n,n3a)
            {24,  {-2, -5, {zai::ALPHA}}},                        // (n,2na)
            {25,  {-2, -6, {zai::ALPHA}}},                        // (n,3na)
            {28,  {-1, -1, {PROTON}}},                            // (n,np)
            {29,  {-4, -8, {zai::ALPHA, zai::ALPHA}}},            // (n,n2a)
            {32,  {-1, -2, {DEUTERON}}},                          // (n,nd)
            {33,  {-1, -3, {TRITON}}},                            // (n,nt)
            {34,  {-2, -3, {HELION}}},                            // (n,nHe3)
            {37,  {0,  -3, {}}},                                  // (n,4n)
            {41,  {-1, -2, {PROTON}}},                            // (n,2np)
            {42,  {-1, -3, {PROTON}}},                            // (n,3np)
            {44,  {-2, -2, {PROTON, PROTON}}},                    // (n,n2p)
            {45,  {-3, -5, {PROTON, zai::ALPHA}}},                // (n,npa)
            {102, {0,  1,  {}}},                                  // (n,g)
            {103, {-1, 0,  {PROTON}}},                            // (n,p)
            {104, {-1, -1, {DEUTERON}}},                          // (n,d)
            {105, {-1, -2, {TRITON}}},                            // (n,t)
            {106, {-2, -2, {HELION}}},                            // (n,He3)
            {107, {-2, -3, {zai::ALPHA}}},                        // (n,a)
            {108, {-4, -7, {zai::ALPHA, zai::ALPHA}}},            // (n,2a)
            {111, {-2, -1, {PROTON, PROTON}}},                    // (n,2p)
            {112, {-3, -4, {PROTON, zai::ALPHA}}},                // (n,pa)
            {113, {-5, -10, {TRITON, zai::ALPHA, zai::ALPHA}}},   // (n,t2a)
            {114, {-5, -9, {DEUTERON, zai::ALPHA, zai::ALPHA}}},  // (n,d2a)
            {115, {-2, -2, {PROTON, DEUTERON}}},                  // (n,pd)
            {116, {-2, -3, {PROTON, TRITON}}},                    // (n,pt)
            {117, {-3, -5, {DEUTERON, zai::ALPHA}}},              // (n,da)
        };

        bool is_reaction_mt(int mt_number);

        ReactionProductT compute_reaction_products(int parent_zai, int mt_number);

        // collapse a histogram (INT=1) cross section in barns with the flux spectrum
        // returns the reaction rate per atom in 1/s
        double collapse_cross_section(const core::Spectrum<double>& in_spectrum,
                                      const std::vector<double>& xs_energies,
                                      const std::vector<double>& xs_values);

        // adds the rate for the parent going to the products
        // uses the same layout as the decay matrix: (parent, daughter)
        void add_reaction_rate(Eigen::MatrixXd& reaction_mat, const std::map<int, int>& zai_map,
                               int parent_zai, int mt_number, double rate);

        // one group-wise (INT=1) cross section in barns, for one parent and one MT
        struct CrossSection {
            int parent_zai;
            int mt_number;
            std::vector<double> energies;
            std::vector<double> values;
        };

        // all the cross sections we can use from a tape
        // read this once, then collapse it with as many spectra as needed
        using ReactionData = std::vector<CrossSection>;

        ReactionData read_reaction_data_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        // convience function to avoid calling ENDFtk fromFile
        ReactionData read_reaction_data_from_file(std::string filename);

        Eigen::MatrixXd construct_reaction_matrix(const ReactionData& reaction_data,
                                                  const core::Spectrum<double>& in_spectrum);

        Eigen::MatrixXd construct_reaction_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                            const core::Spectrum<double>& in_spectrum);
        // convience function to avoid calling ENDFtk fromFile
        Eigen::MatrixXd construct_reaction_matrix_from_file(std::string filename,
                                                            const core::Spectrum<double>& in_spectrum);
    }
}

#endif // ENDF_REACTION_HPP
