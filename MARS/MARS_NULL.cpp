#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<numeric>
#include<boost/program_options>

namespace po = boost::program_options;

template<class T><class F>
std::vector<int> top50(const std::vector<T> &fS, const std::vector<F> fR, const int& topNum)
{
    // initialize original index locations
    std::vector<int> idx(fs.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

    return idx;
}

int main(int argc, char* argv[]){
    boost::program_options::options_description desc("Allowed Options");

    desc.add_options()
    ("help", "Print Help")
    ("n", po::value<int>(), "number of simulations")
    ("g", po::value<std::string>(), "genotype")
    ("o", po::value<std::string>(), "output")
    ("f", po::value<std::string>(), "MARS/fastMARS(0/1, default:0)")
    ("t", po::value<int>(), "topNum(default:50)")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    int fast = 0;
    int topNum = 50;
    int simNum;
    std::string genotypePath, outputPath, fast;
    if(vm.count("help")){
        std::cout << desc <<std::endl;
        exit(1);
    }

    if(vm.count("n")){
        simNum = vm["n"].as<int>();
    }else{
        std::cout << "Please Specify Number of Simulation(n)"<<std::endl;
        std::cout << desc<<std::endl;
        exit(1);
    }
    if(vm.count("g")){
        genotypePath = vm["g"].as<std::string>();
    }else{
        std::cout << "Please Specify genotype path(g)"<<std::endl;
        std::cout << desc << std::endl;
        exit(1);
    }

    if(vm.count("o")){
        outputPath = vm["o"].as<std::string>();
    }else{
        std::cout << "Please Specify output(o)" << std::endl;
        exit(1);
    }

    if(vm.count("f")){
        if(vm["f"].as<std::string>() == "MARS")
            f = 0;
        else if(vm["f"].as<std::string>() == "fastMARS")
            f = 1;
        else{
            std::cout << "Unidentified Option(-f) : " << vm["f"].as<std::vector>()<<std::endl;
            exit(1);
        }
    }else std::cout << "Default MARS/fastMARS(0/1, default:0) is set to 0"<<std::endl;

    if(vm.count("t")){
        topNum = vm["t"].as<int>();
    }else{
        std::cout << "Top Number(t) is set to default(topNum:50)"<<std::endl;
    }

    float** R = (float**)malloc(sizeof(float*)*simNum);
    for(int i=0;i<simNum;i++){
        R[i] = (float*)malloc(sizeof(float)*topNum*2);
        for(int j =0;j<topNum*2;j++)
            R[i][j] = 0.0f;
    }


    time_t start_time = time(NULL);
    






    return 0;
}