using namespace std;

void readData(
                vector<double> &time,
                vector<vector<double>> &pos_x,
                vector<vector<double>> &pos_y,
                vector<vector<double>> &pos_z,
                vector<vector<double>> &d,
                vector<vector<double>> &vol,
                vector<vector<double>> &vel_x,
                vector<vector<double>> &vel_y,
                vector<vector<double>> &vel_z,
                const string &logFile
             )
{




    ifstream Logfile;
    Logfile.open(logFile);
    if(!Logfile.is_open())
      {
          assert(Logfile.is_open());
      }
    string line;
    string letter;
    int ii = 0; //number of data points
    bool t = false;
    while(Logfile.good())
    //for(auto i = 0; i != 10; ++i)
    {
        vector<double> tmpPos_x;
        vector<double> tmpPos_y;
        vector<double> tmpPos_z;

        vector<double> tmpD;
        vector<double> tmpVol;

        vector<double> tmpVel_x;
        vector<double> tmpVel_y;
        vector<double> tmpVel_z;

        unsigned int N = 0;
        getline(Logfile, line);
        stringstream tempLine(line);
        getline(tempLine, letter, ' ');
        //Read current Time
        if (letter == "Time")
        {
            getline(tempLine, letter, ' ');
            getline(tempLine, letter, ' ');
            //cout << letter << endl;
            time.push_back(atof(letter.c_str()));
            ++ii;
        }
        else continue;
        //Read number of dropplets
        getline(Logfile, line);
        stringstream tempLine1(line);
        getline(tempLine1, letter, '\t');
        getline(tempLine1, letter, ';');
        N = atoi(letter.c_str());


        //Read positions
        for(auto j = 0; j != N; ++j)
        {
            getline(Logfile, line);
            stringstream tempLine2(line);
            getline(tempLine2, letter, '(');
            getline(tempLine2, letter, ' ');
            tmpPos_x.push_back(atof(letter.c_str()));
            getline(tempLine2, letter, ' ');
            tmpPos_y.push_back(atof(letter.c_str()));
            getline(tempLine2, letter, ')');
            tmpPos_z.push_back(atof(letter.c_str()));
        }
        //Read Diameter
        getline(Logfile, line);
        for(auto j = 0; j != N; ++j)
        {
            getline(Logfile, line);
            stringstream tempLine2(line);
            getline(tempLine2, letter, '\t');
            getline(tempLine2, letter, ';');
            //cout << letter << endl;
            tmpD.push_back(atof(letter.c_str()));
        }

        //Read Volume
        getline(Logfile, line);
        for(auto j = 0; j != N; ++j)
        {
            getline(Logfile, line);
            stringstream tempLine2(line);
            getline(tempLine2, letter, '\t');
            getline(tempLine2, letter, ';');
            //cout << letter << endl;
            tmpVol.push_back(atof(letter.c_str()));
        }

        //Read velocity
        getline(Logfile, line);
        for(auto j = 0; j != N; ++j)
        {
            getline(Logfile, line);
            stringstream tempLine2(line);
            getline(tempLine2, letter, '(');
            getline(tempLine2, letter, ' ');
            tmpVel_x.push_back(atof(letter.c_str()));
            getline(tempLine2, letter, ' ');
            tmpVel_y.push_back(atof(letter.c_str()));
            getline(tempLine2, letter, ')');
            tmpVel_z.push_back(atof(letter.c_str()));
        }
    cout << "Reading data, iteration - "<<ii <<";" <<endl;
    getline(Logfile, line);

    pos_x.push_back(tmpPos_x);
    pos_y.push_back(tmpPos_y);
    pos_z.push_back(tmpPos_z);
    d.push_back(tmpD);
    vol.push_back(tmpVol);
    vel_x.push_back(tmpVel_x);
    vel_y.push_back(tmpVel_y);
    vel_z.push_back(tmpVel_z);

    }
}

void passedDroplets
        (
            vector<double> &time,
            vector<vector<double>> &pos_x,
            vector<vector<double>> &vel_x,
            vector<vector<double>> &d,
            vector<vector<unsigned int>> &list,
            double surf,
            double delta,
            double &max_d,
            double &d_t
        )
{

    for(auto i = 0; i < time.size(); ++i)
    {
        for(auto j = 0; j < pos_x[i].size(); ++j)
        {
            vector<unsigned int> tmp;
            if ((pos_x[i][j]+(vel_x[i][j]*d_t)) > surf && pos_x[i][j] < surf && d[i][j] < max_d)
            {
                tmp.push_back(i);
                tmp.push_back(j);
                list.push_back(tmp);
            }
            /*
            if (pos_x[i][j] < (surf+delta) && pos_x[i][j] > (surf-delta) && d[i][j] < max_d)
            {
                tmp.push_back(i);
                tmp.push_back(j);
                list.push_back(tmp);
            }
            */
        }
    }
}

vector<vector<double>> radius_zy
        (
            const vector<vector<double>> &pos_y,
            const vector<vector<double>> &pos_z
        )
{
    vector<vector<double>> radius;


    for(auto i = 0; i < pos_y.size(); ++i)
    {
        vector<double> tmpRad;
        for(auto j = 0; j < pos_y[i].size(); ++j)
        {

            double tmp = sqrt((pos_y[i][j]*pos_y[i][j]) +(pos_z[i][j]*pos_z[i][j]));
            tmpRad.push_back(tmp);
        }
        radius.push_back(tmpRad);
    }
    return radius;

}


vector<unsigned int> diameterDistrubution
        (
            const vector<vector<double>> &d,
            const vector<vector<unsigned int>> &list,
            const unsigned int num_div,
            double &d_max,
            double &d_min
        )
{
    unsigned int N[num_div];
    vector<double> dtmp;
    double deltaD;
    vector<unsigned int> diameterDistr;


    //Initialize as zero
    for(unsigned int i = 0; i < num_div; ++i)
    {
        N[i] = 0;
    }
    //create additional arays
    for(auto i = 0; i < list.size(); ++i)
    {
        int j = list[i][0];
        int k = list[i][1];
        dtmp.push_back(d[j][k]);
    }
    //find maximum diameter (assume that min = 0)
    vector<double>::iterator itr_max;
    itr_max = max_element(dtmp.begin(),dtmp.end());
    d_max = dtmp[distance(dtmp.begin(),itr_max)];

    vector<double>::iterator itr_min;
    itr_min = min_element(dtmp.begin(),dtmp.end());
    d_min = dtmp[distance(dtmp.begin(),itr_min)];

    deltaD = (d_max-d_min)/num_div;

    //fill out diameters range
    for(auto i = 0; i < dtmp.size(); ++i)
    {
        int j = (int)((dtmp[i]-d_min)/deltaD);
        N[j] =N[j] + 1 ;
    }

    for(unsigned int i = 0; i < num_div; ++i)
    {
        diameterDistr.push_back(N[i]);
    }

    return diameterDistr;
}
vector<double> diameterCDF
    (
        vector<unsigned int> dDistr
    )
{
    vector<double> diameterCDF;
    unsigned int n = 0;
    for (auto i = 0; i < dDistr.size(); ++i)
    {
        n += dDistr[i];
    }
    for (auto i = 0; i < dDistr.size(); ++i)
    {
        double val = 0.0;
        if(i == 0) val = (double)dDistr[i]/n;//(double)(dDistr[i]/n);
        else val = (double)dDistr[i]/n + diameterCDF[i-1];
        cout << val << endl;
        diameterCDF.push_back(val);
    }
    return diameterCDF;
}
vector<double> scalarAvareging
        (
            const vector<vector<double>> &r,
            const vector<vector<double>> &scalar,
            const vector<vector<unsigned int>> &list,
            const unsigned int &num_div,
            const double &r_max,
            const double &r_min
        )
{

    vector<double> tmpScalar;
    vector<double> tmpR;
    vector<double> res;

    unsigned int N[num_div];
    double S[num_div];
    double deltaR = r_max/num_div;

    //Initialize as zero
    for(unsigned int i = 0; i < num_div; ++i)
    {
        N[i] = 0.0;
        S[i] = 0.0;
    }

    for(auto i = 0; i < list.size(); ++i)
    {
        int j = list[i][0];
        int k = list[i][1];
        tmpScalar.push_back(scalar[j][k]);
        tmpR.push_back(r[j][k]);
    }

    //fill out scalar range
    for(auto i = 0; i < tmpR.size(); ++i)
    {
        int j = (int)((tmpR[i]*0.9999)/deltaR);
        N[j] = N[j] + 1 ;
        S[j] = tmpScalar[i] + S[j];
    }
    //find the avarage
    for(auto i = 0; i < num_div; ++i)
    {
        if (S[i] != 0.0) S[i] = S[i]/N[i];
        res.push_back(S[i]);
    }
    return res;
}

vector<double> positionFromRad
    (
        const double &rad,
        std::mt19937 &generator
    )
{
    double lower_bound = 0;
    double upper_bound = rad;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::uniform_int_distribution<int> dice (0,1);
    int sign1 = dice(generator);
    int sign2 = dice(generator);
    double x = unif(generator);
    double y = sqrt((rad*rad)-(x*x));
    if(sign1 == 0) x = -x;
    if(sign2 == 0) y = -y;
    //cout << "random x - "<< x << endl;
    //cout << "random y - "<< y << endl;
    vector<double> res;
    res.push_back(x);
    res.push_back(y);
    return res;
}
vector<vector<double>> radialVel
    (
        const vector<vector<double>> &pos_y,
        const vector<vector<double>> &pos_z,
        const vector<vector<double>> &vel_y,
        const vector<vector<double>> &vel_z
    )
{

    vector<vector<double>> vrad;
    for (auto i = 0; i < pos_y.size(); ++i)
    {
        vector<double> Vradtmp;
        for (auto j = 0; j < pos_y[i].size(); ++j)
        {
            double alpha = 0;
            double y = pos_y[i][j];
            double z = pos_z[i][j];
            double Vy = vel_y[i][j];
            double Vz = vel_z[i][j];
            if(y > 0.0 && z > 0.0) alpha = atan(abs(z)/abs(y));
            if(y < 0.0 && z > 0.0) alpha = M_PI*0.5 + atan(abs(y)/abs(z));
            if(y < 0.0 && z < 0.0) alpha = M_PI + atan(abs(z)/abs(y));
            if(y > 0.0 && z < 0.0) alpha = M_PI*1.5 +atan(abs(y)/abs(z));

            double v = cos(alpha)*Vy+sin(alpha)*Vz;
            Vradtmp.push_back(v);
        }
        vrad.push_back(Vradtmp);
    }
    return vrad;
}

vector<vector<double>> azimuthalVel
    (
        const vector<vector<double>> &pos_y,
        const vector<vector<double>> &pos_z,
        const vector<vector<double>> &vel_y,
        const vector<vector<double>> &vel_z
    )
{

    vector<vector<double>> vrad;
    for (auto i = 0; i < pos_y.size(); ++i)
    {
        vector<double> Vradtmp;
        for (auto j = 0; j < pos_y[i].size(); ++j)
        {
            double alpha = 0;
            double y = pos_y[i][j];
            double z = pos_z[i][j];
            double Vy = vel_y[i][j];
            double Vz = vel_z[i][j];
            if(y > 0.0 && z > 0.0) alpha = atan(abs(z)/abs(y));
            if(y < 0.0 && z > 0.0) alpha = M_PI*0.5 + atan(abs(y)/abs(z));
            if(y < 0.0 && z < 0.0) alpha = M_PI + atan(abs(z)/abs(y));
            if(y > 0.0 && z < 0.0) alpha = M_PI*1.5 +atan(abs(y)/abs(z));

            double v = cos(alpha)*Vz-sin(alpha)*Vy;
            Vradtmp.push_back(v);
        }
        vrad.push_back(Vradtmp);
    }
    return vrad;
}

void distrWrite
    (
        const vector<double> dDistr,
        const double &d_max,
        const double &d_min,
        const string &dictName
    )
{
    string realDicrName = "resultData/" + dictName;
    std::ofstream file(realDicrName);

    //Head of the dict
    file << "FoamFile" << std::endl << "{"<< std::endl;
    file << '\t' << "version 2.0;" << std::endl;
    file << '\t' << "format ascii;" << std::endl;
    file << '\t' << "class dictionary;" << std::endl;
    file << '\t' << "location \"\";" << std::endl; //change location if needed
    file << '\t' << "object " << dictName << ";" << std::endl; // use your name
    file << "}" << std::endl;

    file << "Div" << '\t' << dDistr.size() << ";" << std::endl;
    file << "D_max" << '\t' << d_max << ";" << std::endl;
    file << "D_min" << '\t' << d_min << ";" << std::endl;

    for (auto j = 0; j <dDistr.size(); ++j)
    {
         file << "Distr" << j << '\t'<< dDistr[j] << ";" << endl;
    }
}
void gnuWrite
    (
        const vector<vector<double>> &x,
        const vector<vector<double>> &y,
        const vector<vector<unsigned int>> &list,
        const string &name
    )
{
    string realName = "gnuplotData/" + name;
    std::ofstream file(realName);
    for (auto i = 0; i <list.size(); ++i)
    {
         file << x[list[i][0]][list[i][1]] << '\t' << y[list[i][0]][list[i][1]] << endl;
    }
}
void gnuWrite
    (
        const vector<double> &x,
        const vector<double> &y,
        const string &name
    )
{
    string realName = "gnuplotData/" + name;
    std::ofstream file(realName);
    for (auto i = 0; i <x.size(); ++i)
    {
         file << x[i] << '\t' << y[i] << endl;
    }
}

vector<double> readPostData
    (
        const string &fileName
    )
{
    vector<double> data;
    std::fstream file;
    file.open(fileName);
    if(!file.is_open())
      {
          assert(file.is_open());
      }
    std::string line;
    std::string letter;
    while(file.good())
    //for (int var = 0; var < 10; ++var)
    {
        getline(file, line);
        std::stringstream tempLine(line);
        getline(tempLine, letter, ';');
        //std::cout << letter  << endl;
        if(letter != "") data.push_back(std::stof(letter));
        else break;
    }
    return data;
}

void gnuWritePDF
    (
        const string &fileName,
        const vector<double> &scalar,
        const int Div
    )
{
    double min = 1e5;
    double max = -1e5;
    double tmp[Div];

    for (auto i = 0; i < Div; ++i)
    {
        tmp[i] = 0.0;
    }

    for (auto i = 0; i < scalar.size(); ++i)
    {
        if(min > scalar[i]) min = scalar[i];
        if(max < scalar[i]) max = scalar[i];
    }
    for (auto i = 0; i < scalar.size(); ++i)
    {
        int j = ((scalar[i] - min)/(max-min))*Div;
        tmp[j] += 1;
    }

    std::ofstream file(fileName);
    for (auto i = 0; i < Div; ++i)
    {
         file << i << '\t' << tmp[i] << endl;
    }
}

double meanValuePDF (const vector<unsigned int>& pdf,const double &d_max)
{
    double res = 0.0;
    double tmp = 0.0;
    for (auto i = 0; i < pdf.size(); ++i)
    {
        res += (double)pdf[i];
    }
    res = res / 2;
    for (auto i = 0; i < pdf.size(); ++i)
    {
        tmp += (double)pdf[i];
        if(tmp > res) return (((double)i/(double)pdf.size())*d_max);
    }
    return 0.0;
}

float momentumPlane
(
    vector<vector<double>> &vel_x,
    vector<vector<double>> &vol,
    vector<vector<unsigned int>> &list,
    double &rho
)
{
    int ii_old = 0;
    vector <double> tmpVec;
    double tmp = 0;
    for (auto i = 0; i < list.size(); ++i)
    {
        int ii = list[i][0];
        int jj = list[i][1];
        //cout << ii << endl;
        if(ii_old == ii)
        {

            tmp += vel_x[ii][jj]*vol[ii][jj]*rho;
        }
        else
        {
            tmpVec.push_back(tmp);

            tmp = vel_x[ii][jj]*vol[ii][jj]*rho;
            ++ii_old;
        }
    }
    float res = std::accumulate(tmpVec.begin(), tmpVec.end(), 0.0)/tmpVec.size();
    return res;

}
