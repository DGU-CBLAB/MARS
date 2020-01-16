#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<numeric>
#include <RcppArmadillo.h>
#include <mvnorm.h>
#include<boost/program_options>

namespace po = boost::program_options;
template<class T>
void transpose(vector<vector<T> > &b)
{
    if (b.size() == 0)
        return;

    vector<vector<T> > trans_vec(b[0].size(), vector<T>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    b = trans_vec;    // <--- reassign here
}

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



    time_t start_time = time(NULL);

    // Variables
    float** R;
    std::vector<std::vector<std::string>>> x;
    int m, n, d;

    // R matrix
    R = (float**)malloc(sizeof(float*)*simNum);
    for(int i=0;i<simNum;i++){
        R[i] = (float*)malloc(sizeof(float)*topNum*2);
        for(int j =0;j<topNum*2;j++)
            R[i][j] = 0.0f;
    }

    // x matrix(genotype file - as whole)
    std::string str;
    try{
        ofstream genotype_file(genotypePath);
        while(std::readLine(genotype_file, str)){
            std::vector<std::string> vec_str = std::split(str,'\t');
            x.push_back(vec_str);
        }
    }catch(std::Exception& e){
        std::cout<< e.what() << std::endl;
        exit(1);
    }

    // m, n
    m = x.size();
    n = x.at(0).size();

    // x = x[,2:n];
    for(int i=0; i<x.size(); i++){
        vector<std::string>::const_iterator first = x.at(i).begin() + 3;
        vector<std::string>::const_iterator last = x.at(i).end();
        std::vector<std::string> replace(first, last);
        x.at(i).clear();
        x.at(i) = replace;
    }

    n = n-1;
    d = x.size();

    // x <- as.numeric(x)
    vector<vector<float>> x_;
    for(int i =0;i<x.size();i++){
        vector<float> d_v;
        for(int j =0;j<x.at(i).size();j++){
            d_v.push_back(stof(x.at(i).at(j)));
        }
        x_.push_back(d_v);
    }

    transpose(x_);

    // xs = scale(x,center=TRUE, scale=TRUE)
    vector<vector<float>> xs;

    vector<float> v;
    for(int i=0;i<x_.size();i++){
        for(int j=0;j<x_.at(i).size();j++){
            total.push_back(x_.at(i).(j));
        }
    }
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);


    for(int i=0;i<x_.size();i++){
        vector<float> tmp;
        for(int j=0;j<x_.at(i).size();j++){
            tmp.push_back(x_.at(i).at(j)-stdev);
        }
        xs.push_back(tmp);
    }
    float avg = sum / cnt;

    // I = diag(n)
    vector<vector<float>> I;
    for(int i =0;i<n;i++){
        vector<float> tmp;
        for(int j=0;j<n;j++){
            if(i==j)
                tmp.push_back(i+1);
            else
                tmp.push_back(0);
        }
        I.push_back(tmp);
    }

    // mu = mat.or.vec(nr, nc)
    vector<float> mu;
    for(int i=0;i<n;i++)
        mu.push_back(0);








    if(fast == 0){
        std::cout << "Generating null for MARS" << std::endl;
        Sall = rmvnorm(simNum, mu, I);

    }else{
        std::cout << "Generating null for fastMARS" << std::endl;





    }


    return 0;
}