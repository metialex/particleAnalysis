using namespace std;

class scalarPDF
{
public:
    scalarPDF
        (
            const vector<vector<double>> &r,
            const vector<vector<double>> &scalar,
            const vector<vector<unsigned int>> &list,
            const int &scalarAvSize
        )
    {
        vector<double> x_data,y_data; //x correspond to R, y to scalar
        //create additional arays
        for(auto i = 0; i < list.size(); ++i)
        {
            int j = list[i][0];
            int k = list[i][1];
            x_data.push_back(r[j][k]);
            y_data.push_back(scalar[j][k]);
        }
        //find min/max R
        min_r = 1e5;
        max_r = -1e5;
        for(auto i = 0; i < x_data.size(); ++i)
        {
            if (min_r > x_data[i]) min_r = x_data[i];
            if (max_r < x_data[i]) max_r = x_data[i];
        }
        double range = max_r - min_r;
        //fill out data arrays

        for(auto i = 0; i < scalarAvSize; ++i)
        {
            vector<double> tmpData_x;
            vector<double> tmpData_y;
            for(auto j = 0; j < x_data.size(); ++j)
            {
                if(x_data[j] > (min_r + i*(range/scalarAvSize)) && x_data[j] < (min_r + (i+1)*(range/scalarAvSize)))
                {
                    tmpData_x.push_back(x_data[j]);
                    tmpData_y.push_back(y_data[j]);
                }
            }
            this->data_x.push_back(tmpData_x);
            this->data_y.push_back(tmpData_y);
        }
        //fill out min/max vector
        for(auto i = 0; i < data_y.size(); ++i)
        {
            double min = 1e5;
            double max = -1e5;
            vector<double> tmpMinMax;
            for(auto j = 0; j < data_y[i].size(); ++j)
            {
                if (min > data_y[i][j]) min = data_y[i][j];
                if (max < data_y[i][j]) max = data_y[i][j];
            }
            tmpMinMax.push_back(min);
            tmpMinMax.push_back(max);
            minMax.push_back(tmpMinMax);
        }


    }

    scalarPDF
        (
            const vector<vector<double>> &r,
            const vector<vector<double>> &scalar,
            const vector<vector<unsigned int>> &list,
            const int &scalarAvSize,
            const double &min
        )
    {
        vector<double> x_data,y_data; //x correspond to R, y to scalar
        //create additional arays
        for(auto i = 0; i < list.size(); ++i)
        {
            int j = list[i][0];
            int k = list[i][1];
            x_data.push_back(r[j][k]);
            y_data.push_back(scalar[j][k]);
        }
        //find min/max R
        min_r = min;
        max_r = -1e5;
        for(auto i = 0; i < x_data.size(); ++i)
        {
            if (max_r < x_data[i]) max_r = x_data[i];
        }
        double range = max_r - min_r;
        //fill out data arrays

        for(auto i = 0; i < scalarAvSize; ++i)
        {
            vector<double> tmpData_x;
            vector<double> tmpData_y;
            for(auto j = 0; j < x_data.size(); ++j)
            {
                if(x_data[j] > (min_r + i*(range/scalarAvSize)) && x_data[j] < (min_r + (i+1)*(range/scalarAvSize)))
                {
                    tmpData_x.push_back(x_data[j]);
                    tmpData_y.push_back(y_data[j]);
                }
            }
            this->data_x.push_back(tmpData_x);
            this->data_y.push_back(tmpData_y);
        }
        //fill out min/max vector
        for(auto i = 0; i < data_y.size(); ++i)
        {
            double min = 1e5;
            double max = -1e5;
            vector<double> tmpMinMax;
            for(auto j = 0; j < data_y[i].size(); ++j)
            {
                if (min > data_y[i][j]) min = data_y[i][j];
                if (max < data_y[i][j]) max = data_y[i][j];
            }
            tmpMinMax.push_back(min);
            tmpMinMax.push_back(max);
            minMax.push_back(tmpMinMax);
        }


    }

    void createPDF(const unsigned int &num_div)
    {
        //create pdfSum vector
        for(auto i = 0; i < data_y.size(); ++i)
        {
            pdfSum.push_back(data_y[i].size());
        }
        //create pdf vector
        for(auto i = 0; i < data_y.size(); ++i)
        {
            double range = minMax[i][1] - minMax[i][0];
            unsigned int N[num_div];
            for(auto j = 0; j < num_div; ++j) N[j] = 0;
            vector<unsigned int> tmpPdf;
            for(auto j = 0; j < data_y[i].size(); ++j)
            {
                for(auto k = 0; k < num_div; ++k)
                {
                    if(data_y[i][j] >= minMax[i][0] + k*(range/num_div) && data_y[i][j] <= 1.00001*(minMax[i][0] + (k+1)*(range/num_div))) N[k] += 1;  //minMax[i][0] - min, [i][1] - max
                }
            }
            for(auto j = 0; j < num_div; ++j)
            {
                tmpPdf.push_back(N[j]);
            }
            pdf.push_back(tmpPdf);
        }
    }
    void createCDF()
    {
        for(auto i = 0; i < this->pdfSum.size(); ++i)
        {
            vector<double> tmpCdf;
            for(auto j = 0; j < this->pdf[i].size(); ++j)
            {
                double x = 0;
                if(j == 0) x = (double)(this->pdf[i][j])/(double)(this->pdfSum[i]);
                else x = (double)(this->pdf[i][j])/(double)(this->pdfSum[i]) + tmpCdf[j-1];
                tmpCdf.push_back(x);
            }
            this->cdf.push_back(tmpCdf);
        }
    }
    void printData()
    {
        cout << "data size - " << this->data_x.size() << endl << endl;
        //print X values
        for(auto i = 0; i < this->data_x.size(); ++i)
        {
            cout << i << " - "<< this->data_x[i].size() <<endl;
            for(auto j = 0; j < this->data_x[i].size(); ++j)
            {
                cout << this->data_x[i][j] <<endl;
            }
        }

        //print Y values
        for(auto i = 0; i < this->data_y.size(); ++i)
        {
            cout << i << " - "<< this->data_y[i].size() <<endl;
            for(auto j = 0; j < this->data_y[i].size(); ++j)
            {
                cout << this->data_y[i][j] <<endl;
            }
        }
        /*
        //print min max
        for(auto i = 0; i < this->minMax.size(); ++i)
        {
            cout << i << endl;
            for(auto j = 0; j < this->minMax[i].size(); ++j)
            {
                cout << this->minMax[i][j] <<endl;
            }
        }
        */
        //print pdf
        for(auto i = 0; i < this->pdf.size(); ++i)
        {
            cout << i << " - "<< this->pdf[i].size() <<endl;
            for(auto j = 0; j < this->pdf[i].size(); ++j)
            {
                cout << this->pdf[i][j] <<endl;
            }
        }
        //print pdfSum
        for(auto i = 0; i < this->pdfSum.size(); ++i)
        {
            cout << i << " - " << pdfSum[i] << endl;
        }
        //print cdf
        for(auto i = 0; i < this->cdf.size(); ++i)
        {
            cout << i << " - "<< this->cdf[i].size() <<endl;
            for(auto j = 0; j < this->cdf[i].size(); ++j)
            {
                cout << this->cdf[i][j] <<endl;
            }
        }
    }
    void createGlobalPDF
        (
            const vector<vector<double>> &scalar,
            const vector<vector<unsigned int>> &list,
            const int &Div
        )
    {
        //find min/max of the range
        double min = 1e5;
        double max = -1e5;
        for (auto i = 0; i < list.size(); ++i)
        {
            if (min > scalar[list[i][0]][list[i][1]]) min = scalar[list[i][0]][list[i][1]];
            if (max < scalar[list[i][0]][list[i][1]]) max = scalar[list[i][0]][list[i][1]];
        }
        //distrubute according to the velosity
        double tmp[Div];
        for (auto i = 0; i < Div; ++i)
        {
            tmp[i] = 0;
        }
        for (auto i = 0; i < list.size(); ++i)
        {
            int j = (((scalar[list[i][0]][list[i][1]]) - min)/(max - min))*Div;
            tmp[j] += 1;
        }
        for (auto i = 0; i < Div; ++i)
        {
            this->globalPDF.push_back(tmp[i]);
        }
    }
    void saveGlodablPDF(const string &dictName)
    {
        string realDictName = "resultData/" + dictName;
        std::ofstream file(realDictName);
        for (auto i = 0; i <globalPDF.size(); ++i)
        {
             file << i << '\t' << globalPDF[i] << endl;
        }
    }
    void saveCDF(const string &dictName)
    {
        string realDictName = "resultData/" + dictName;
        std::ofstream file(realDictName);

        //Head of the dict
        file << "FoamFile" << std::endl << "{"<< std::endl;
        file << '\t' << "version 2.0;" << std::endl;
        file << '\t' << "format ascii;" << std::endl;
        file << '\t' << "class dictionary;" << std::endl;
        file << '\t' << "location \"\";" << std::endl; //change location if needed
        file << '\t' << "object " << dictName << ";" << std::endl; // use your name
        file << "}" << std::endl;

        //number of divisions
        file << "Div" << '\t' << cdf.size() << ";" << std::endl;
        //number of subdivisions
        file << "subDiv" << '\t' << cdf[0].size() << ";" << std::endl << std::endl;
        //min max R
        file << "min_r" << '\t' << min_r << ";" <<endl;
        file << "max_r" << '\t' << max_r << ";" <<endl;

        //print the cdf itself
        for (auto i = 0; i < cdf.size(); ++i)
        {
            file << "localMin"<< i << '\t' << minMax[i][0]<< ";" << std::endl;
            file << "localMax"<< i << '\t' << minMax[i][1]<< ";" << std::endl;
            //file << "(";
            for (auto j = 0; j < cdf[i].size(); ++j)
            {
                 file<< "cdf"<< i<< "_" << j<< '\t' << cdf[i][j] << ";" << endl;
            }
            //file << ")"<<std::endl;
        }
    }

    void savePDF(const string &dictName)
    {
        //print the pdf itself
        for (auto i = 0; i < cdf.size(); ++i)
        {
            string realDictName = "resultData/" + dictName + std::to_string(i);
            std::ofstream file(realDictName);
            for (auto j = 0; j < cdf[i].size(); ++j)
            {
                //refers to radius values
                double ii = min_r + ((double)i/cdf.size())*(max_r-min_r);
                //refers to scalar values
                double jj = minMax[i][0] + (((double)j/cdf[i].size()))*(minMax[i][1] - minMax[i][0]);
                 file << ii << " " << jj << " " << pdf[i][j]  << endl;
            }
        }
    }

    void saveCDFNew(const string &dictName)
    {
        //print the pdf itself
        for (auto i = 0; i < cdf.size(); ++i)
        {
            string realDictName = "resultData/" + dictName + std::to_string(i);
            std::ofstream file(realDictName);
            double iii = min_r + ((double)i/cdf.size())*(max_r-min_r);
            file << iii << " 0 0" << endl;
            for (auto j = 0; j < cdf[i].size(); ++j)
            {
                //refers to radius values
                double ii = min_r + ((double)i/cdf.size())*(max_r-min_r);
                //refers to scalar values
                double jj = minMax[i][0] + (((double)j/cdf[i].size()))*(minMax[i][1] - minMax[i][0]);
                 file << ii << " " << j+1 << " " << cdf[i][j]  << endl;
            }
        }
    }
    ~scalarPDF(){}
private:
    vector<vector<double>> data_x, data_y, minMax;
    vector<vector<unsigned int>> pdf;
    vector<vector<double>> cdf;
    vector<unsigned int> pdfSum;
    double min_r, max_r;
    vector<double> globalPDF;
};
