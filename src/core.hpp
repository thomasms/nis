#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "nuclides.hpp"
#include "util.hpp"

namespace nis {
    namespace core {

        const double LN_2 = std::log(2);

        double compute_decay_const(double half_life){
            return LN_2/half_life;
        }

        // class IProductComputer{
        //         std::vector<int> compute() const = 0;
        //     private:
        //         int _rtyp;
        //         int _parent_zai;
        // };

        // class GammaRayReactionProducts : public IProductComputer{
        //     GammaRayReactionProducts(int parent_zai) : _rtyp(0), _parent_zai(parent_zai) {};
            
        //     std::vector<int> compute() const override{
        //         std::vector<int> products = {_parent_zai};
        //         return products;
        //     };

        // };

        // class BetaDecayReactionProducts : public IProductComputer{
        //     BetaDecayReactionProducts(int parent_zai) : _rtyp(1), _parent_zai(parent_zai) {};
            
        //     std::vector<int> compute() const override{
        //         std::vector<int> products = {_parent_zai+10, };
        //         return products;
        //     };

        // };

        // std::vector<int> compute_products(int parent_zai, int rtyp){
        // }


        /*
        Needs to be cleaned up and move somewhere else.
        */
        void _update_decay_matrix(const njoy::ENDFtk::tree::Material& material, const std::map<int, int>& zai_map, Eigen::MatrixXd& decay_mat){
            constexpr int file_number = 8;
            constexpr int mt_number = 457;
            if(material.hasMF(file_number)){
                njoy::ENDFtk::tree::File file = material.file(file_number);
                if(file.hasMT(mt_number)){
                    njoy::ENDFtk::tree::Section decay_data = file.MT(mt_number);
                    auto parsed_section = decay_data.parse<file_number, mt_number>();
                    const int ZA = parsed_section.ZA();
                    const int ZAI = ZA*10;
                    if (zai_map.count(ZAI) == 1){
                        const int za_index = zai_map.at(ZAI);
        
                        if( !parsed_section.isStable() ){
                            const double half_life = parsed_section.halfLife()[0];
                            const double decay_const = compute_decay_const(half_life);

                            // the diagonal is the rate lost from the parent
                            decay_mat(za_index, za_index) -= decay_const;

                            // todo: loop over decay states and compute the outgoing
                            // for(auto decay_mode: parsed_section.decayModes()){
                                
                            // }
                        }
                    }
                }
            }
        }

        Eigen::MatrixXd construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape){

            // this could be done at compile time using std::array
            const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);
            // add 1 for the sink
            const int mat_nr_rows = zai_map.size() + 1;
            const int mat_nr_cols = mat_nr_rows;

            Eigen::MatrixXd decay_mat = Eigen::MatrixXd::Zero(mat_nr_rows, mat_nr_cols);

            const std::vector<int> mat_numbers = tape.materialNumbers();
            // loop over all materials in the file, most likely just one
            for(auto mat_number: mat_numbers){
                // this can return multiple instances of the same MAT number
                // due to different temperatures
                // loop over temperatures for each MAT
                auto materials = tape.material(mat_number);
                for(auto material: materials){
                    _update_decay_matrix(material, zai_map, decay_mat);
                }
            }

            return decay_mat;
        }


    }
}