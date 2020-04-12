    //Variables for initial data
    vector<double> time;
    vector<vector<double>> pos_x;
    vector<vector<double>> pos_y;
    vector<vector<double>> pos_z;
    vector<vector<double>> d;
    vector<vector<double>> vol;
    vector<vector<double>> vel_x;
    vector<vector<double>> vel_y;
    vector<vector<double>> vel_z;
    vector<vector<double>> vel_Rad;
    vector<vector<double>> vel_Az;
    vector<vector<double>> r;
    vector<vector<double>> momentum_x;
    string logFile = "initialData/output0.5-0.45";
    //variable for random gen al
    std::mt19937 generator;
    generator.seed(std::time(0));

    //Variables for analisys
    vector<vector<unsigned int>> list; //first value list[i][0] - corresponds to the index (time),
                              //the second to index of dropplet list[i][1]
    double surf =  0.003; // value correspond to the x coordinate of the control plane
    double delta = 0.00005; //correspond to the max distance from CP where droplets are considered
    const double r_min = 0.0;
    const double r_max = 0.0004;
    unsigned int num_div = 20;
    int pdfSize = 20;
    double d_max = 0; //changing during diameterDistrubution
    double d_min = 0; //changing during diameterDistrubution
    double dd_max = 2e-5;
    double rho = 148.0;
    double d_t = 1e-7;
