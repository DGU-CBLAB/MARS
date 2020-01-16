#include<iostream>
#include<fstream>
#include<numeric>
#include<vector>
#include<algorithm>
#include<boost/program_options.hpp>

namespace po = boost::program_options;
template <typename T>
std::vector<int> sort_indexes(const std::vector<T> &v)
{

    // initialize original index locations
    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}
int main(int argc, char* argv[]){
    boost::program_options::options_description desc("Allowed Options");
    desc.add_options()
    ("help", "Print Help")
    ("g", po::value<std::string>(), "genotypePath")
    ("s", po::value<std::string>(), "statPath")
    ("o", po::value<std::string>(), "output_genotype")
    ("u", po::value<std::string>(), "output_stat")
    ("t", po::value<int>(),"topNum(default:50)")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    std::string genotypePath, statPath, output_genotype, output_stat;
    int topNum = 50;

    if(vm.count("help")){
        std::cout << desc << std::endl;
        exit(1);
    }

    if(vm.count("g")) genotypePath = vm["g"].as<std::string>();
    if(vm.count("s")) statPath = vm["s"].as<std::string>();
    if(vm.count("o")) output_genotype = vm["o"].as<std::string>();
    if(vm.count("u")) output_stat = vm["u"].as<std::string>();

    try{
        std::ifstream stat_file(statPath);
        std::ifstream geno_file(genotypePath);
        std::string str;
        std::vector<double> S;
        std::vector<std::string> x;
        while(std::getline(stat_file,str)){
            S.push_back(abs(stoi(str)));
        }stat_file.close();
        while(std::getline(geno_file,str)){
            x.push_back(str+"\n");
        }geno_file.close();

        // sort(S.begin(), S.end(), std::greater<>());
        std::vector<int> indices = sort_indexes<double>(S);
        

        std::ofstream out_stat(output_stat);
        for(int i : indices){
            out_stat.write(std::to_string(S.at(i)).c_str(),sizeof(double));
        }
        out_stat.write("\n",sizeof(std::string));
        out_stat.close();

        std::ofstream out_geno(output_genotype);
        for(int i : indices){
            out_geno.write(x.at(i).c_str(),sizeof(std::string));
        }
        out_geno.write("\n", sizeof(std::string));
        out_geno.close();

    }catch(std::exception& e){
        std::cout << e.what() << std::endl;
        exit(-1);
    }
    return 1;
}