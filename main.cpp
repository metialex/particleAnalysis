#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <algorithm>
#include <random>
#include <ctime>
#include "distributionFunc.h"
#include "myFunctions.h"



using namespace std;

int main(int argc, char *argv[])
{
    #include "createVars.H"

    //read Data 
    readData(time,pos_x,pos_y,pos_z,d,vol,vel_x,vel_y,vel_z,logFile);
    //Make Azimuthal and radial Velosities
    vel_Rad = radialVel(pos_y,pos_z,vel_y,vel_z);
    vel_Az = azimuthalVel(pos_y,pos_z,vel_y,vel_z);
    //Momentum X
    for (auto i = 0; i <vel_x.size(); ++i)
    {
        vector<double> tmpRes;
        for (auto j = 0; j <vel_x[i].size(); ++j)
        {
            double tmp = vel_x[i][j]*vol[i][j]*rho;
            tmpRes.push_back(tmp);
        }
        momentum_x.push_back(tmpRes);
    }


    //make a list
    passedDroplets(time,pos_x,vel_x,d,list,surf,delta, dd_max, d_t);

    float moment = momentumPlane(vel_x, vol, list, rho);
    //make radius vector
    r = radius_zy(pos_y,pos_z);

    //calculate diameter distribution
    vector<unsigned int> dDistr = diameterDistrubution(d,list,num_div, d_max, d_min);
    double mean = meanValuePDF (dDistr,d_max);
    cout << "mean D - " << mean << endl;
    cout << "accumulated moment - " << moment << endl;
    vector<double> dCDF = diameterCDF(dDistr);
    distrWrite(dCDF, d_max, d_min, "diamDistrDict");

    //calculate scalar dependency from radius
    vector<double> rrr;
    for (auto i = 0; i< 50; ++i) {rrr.push_back(r_min + (i*(r_max-r_min)/50));}
    vector<double> scalarAv = scalarAvareging(r,momentum_x,list,50,r_max,r_min);
    gnuWrite(rrr,scalarAv,"momentAv" );

    scalarPDF dPDF(d, r, list, num_div, d_min*1.2); //<-check this out!!!
    scalarPDF vel_xPDF(r, vel_x, list, num_div);
    scalarPDF vel_yPDF(r, vel_y, list, num_div);
    scalarPDF vel_zPDF(r, vel_z, list, num_div);
    scalarPDF vel_RadPDF(r, vel_Rad, list, num_div);
    scalarPDF vel_AzPDF(r, vel_Az, list, num_div);
    scalarPDF exmpl(r, d, list, num_div);


    dPDF.createPDF(pdfSize);
    vel_xPDF.createPDF(pdfSize);
    vel_yPDF.createPDF(pdfSize);
    vel_zPDF.createPDF(pdfSize);
    vel_RadPDF.createPDF(pdfSize);
    vel_AzPDF.createPDF(pdfSize);
    exmpl.createPDF(pdfSize);

    dPDF.createCDF();
    vel_xPDF.createCDF();
    vel_yPDF.createCDF();
    vel_zPDF.createCDF();
    vel_RadPDF.createCDF();
    vel_AzPDF.createCDF();
    exmpl.createCDF();

    dPDF.saveCDF("dCDFDict");
    exmpl.savePDF("exmpl");
    vel_xPDF.savePDF("velXPDFDict");
    vel_RadPDF.savePDF("velRadPDFDict");
    vel_xPDF.saveCDF("velXCDFDict");
    vel_yPDF.saveCDF("velYCDFDict");
    vel_zPDF.saveCDF("velZCDFDict");
    vel_RadPDF.saveCDF("velRadCDFDict");
    vel_AzPDF.saveCDF("velAzCDFDict");

    vel_xPDF.saveCDFNew("velXCDFNew");

    //additional staff
    vel_xPDF.createGlobalPDF(vel_x, list, 50);
    vel_xPDF.saveGlodablPDF("velXg");

    vel_RadPDF.createGlobalPDF(vel_Rad, list, 50);
    vel_RadPDF.saveGlodablPDF("velRadg");

    vel_AzPDF.createGlobalPDF(vel_z, list, 50);
    vel_AzPDF.saveGlodablPDF("velAzg");

    dPDF.createGlobalPDF(d, list, 50);
    dPDF.saveGlodablPDF("dg");

    exmpl.createGlobalPDF(d, list, 50);
    exmpl.saveGlodablPDF("explg");

    gnuWrite(r, d, list, "ddd");
    gnuWrite(r, momentum_x, list, "momentumX");
    gnuWrite(r, vel_Rad, list, "velRad");
    gnuWrite(r, vel_Az, list, "velAz");
    gnuWrite(r, d, list, "example");




    return 0;
}
