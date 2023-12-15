/**
 MIT License
 
 Copyright (c) 2021 CNRS
 
 Loïs Paulin, David Coeurjolly, Jean-Claude Iehl, Nicolas Bonneel, Alexander Keller, Victor Ostromoukhov
 "Cascaded Sobol' Sampling", ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia), 40(6), pp. 274:1–274:13, December 2021
 December 2021
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */
#include <iostream>
#include <string>
#include <fstream>
#include <random>

#include "CLI11.hpp"

#include "Samplers/SobolGenerator1D.h"
#include "Samplers/OwenScrambling.h"


using namespace std;

template< typename T >
inline void getCascadedSobol( T *values,
                const std::vector<SobolGenerator1D >& sobols,
                const uint32_t *seeds,
                const int nDims,
                const uint32_t n,
                const uint32_t nbits,
                const bool owen_permut_flag = true,
                const uint32_t owen_tree_depth = 32 )
 {
    uint32_t index = n;		// dim 0: take n-th point as index -> into 32-bit integer
    for(int idim = 0; idim < nDims; idim++)
    {
        index = sobols[idim].getSobolInt(index);	// radix-inversion + sobol
        uint32_t result = index;
        if(owen_permut_flag)						// we need to apply OwenScrambling only when this flag is set
            result = OwenScrambling(result, seeds[idim], owen_tree_depth);
        values[idim] = T(result) / T(std::numeric_limits<uint32_t>::max());	// final value (double) for this dimension
        index = index >> (32-nbits);				// this will be used as new index for the next dimension
    }
}   // getCascadedSobol

int main(int argc, char **argv) {

    //Call parameters handling
    CLI::App app { "cascadedSobol" };
    int nDims = 4;
    app.add_option("-d,--nDims", nDims, "number of dimensions to generate, default: " + std::to_string(nDims)  );
    int npts = 1024;
    app.add_option("-n,--npts", npts, "number of points to generate, default: " + std::to_string(npts) );
    int seed = 133742;
    app.add_option("-s,--seed", seed, "Random number generator seed. default: " + std::to_string(seed)  );
    string input_dir_vectors = "../data/cascaded_sobol_init_tab.dat";
    app.add_option("-i,--idv", input_dir_vectors, "input sobol initialisation (ascii file), default: " + input_dir_vectors );
    bool owen_permut_flag = false;
    app.add_flag("-p", owen_permut_flag, "Apply Owen permutation on output points, default: " + std::to_string(owen_permut_flag) );
    int nbReal = 1;
    app.add_option("-m", nbReal, "number of realizations of the sampler, default: " + std::to_string(nbReal));
    std::string output_fname = "out.dat";
    app.add_option("-o,--output", output_fname, "output fname, default: " + output_fname );
    CLI11_PARSE(app, argc, argv)

    //Create sobol matrices data
    std::vector<SobolGenerator1D > sobols;	// array of sobol data per dim
    std::vector<uint32_t> d;
    std::vector<uint32_t> s;
    std::vector<uint32_t> a;
    std::vector<std::vector<uint32_t>> m;

    //Read sobol matrices value from file
    std::ifstream tableFile(input_dir_vectors);
    if (!tableFile.is_open()) {
        std::cerr << input_dir_vectors << " cannot be read." << std::endl;
        exit(1);
    };
    load_init_table(tableFile, d, s, a, m, nDims);
    init_sobols(sobols, d, s, a, m);

    //Open output file
    std::ofstream ofs(output_fname , std::ios_base::out);

    //init random generator for owen scrambling
    std::uniform_int_distribution<uint32_t> unif32bits(0);  // uniform distribution of uint32_ts between 0 and 2^32-1
    mt19937_64 gen(seed);

    auto nbits = (uint32_t) ((double) log((double) npts) / (double) log((double) 2));
    vector<double> points(npts * nDims);
    //For each realizations
    for(int realization=0; realization < nbReal; ++realization )
    {
        //init owen scrambling seeds
        vector<uint32_t> realSeeds(nDims);
        for (int iDim = 0; iDim < nDims; ++iDim) {
            realSeeds[iDim] = unif32bits(gen);
        }

        //get each point
        for (int ipt = 0; ipt < npts; ipt++) {
            getCascadedSobol(points.data() + ipt * nDims, sobols, realSeeds.data(), nDims, ipt, nbits, owen_permut_flag);
        }

        //export points in ASCII
        for (int ipt = 0; ipt < npts; ipt++) {
            for (unsigned int idim = 0; idim < nDims; idim++) {
                ofs << std::setprecision(20) << std::fixed << points[ipt * nDims + idim] << " ";  // dims {0, ..., nDims-1}
            }
            ofs << endl;
        }
        //Diffent realizations are separated by #
        if (realization+1 != nbReal){
            ofs<<"#"<<std::endl;
        }
    }
    ofs.close();
}

