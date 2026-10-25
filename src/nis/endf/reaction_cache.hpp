#ifndef ENDF_REACTION_CACHE_HPP
#define ENDF_REACTION_CACHE_HPP

#include <string>

#include "nis/endf/reaction.hpp"

namespace nis {
    namespace endf {

        // A binary copy of the reaction data, so we do not parse the ENDF tape (~1 GB) each run.
        // Each energy grid is written once, and each cross section only from its first to its
        // last non-zero value. The numbers are in the byte order of the machine that wrote it.
        void write_reaction_cache(const ReactionData& reaction_data, const std::string& filename);
        ReactionData read_reaction_cache(const std::string& filename);

        // true if the file starts with the cache header
        bool is_reaction_cache(const std::string& filename);

        // reads a cache or an ENDF tape, from what the file is
        ReactionData read_reaction_data(const std::string& filename);
    }
}

#endif // ENDF_REACTION_CACHE_HPP
