/**
 * @file l051.cpp
 * @author Vivian Feng
 * @brief Lab 5, Part 1 Canny Edge detection with Sobel and single threshold
 * @version 0.1
 * @date 2022-12-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <sstream>
#include <float.h>
#include <list>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <iterator>
#include <stack>

using namespace std::chrono;
using namespace std;

const int PRECISION = 23;
const double TAN225 = 0.41421356237;
const double TAN675 = 2.41421356237;

// typedef std::numeric_limits< double > dbl;
// ofstream results;
class Grid
{
private:
    int width;
    int height;
    vector<vector<int>> grid;

public:
    Grid(int width, int height)
    {
        grid = vector<vector<int>>(height, vector<int>(width, 0));
        //grid = new int *[mySize];
        //size = mySize;
        this->height = height;
        this->width=width;
    }
    Grid(const Grid &gr)
    {
        grid = gr.grid;
        width = gr.width;
        height = gr.height;
    }
    ~Grid()
    {
        // cout << "Grid deleted"<<endl;
    }
    void setPixel(int r, int c, int val){
        grid[r][c] = val;
    }
    int getPixel(int r, int c){
        return grid[r][c];
    }
    int getHeight(){
        return height;
    }
    int getWidth(){
        return width;
    }
    void toPPM(string filename,int maxPix)
    {
        ofstream myfile;
        myfile.open(filename);
        myfile << "P3"<< " " << width << " " << height << " " << maxPix << endl;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                myfile << grid[i][j] << " "<< grid[i][j] << " "<<grid[i][j]<<" ";
            }
            myfile << endl;
        }
        myfile.close();
    }
};

/**
 * Given two points a and b which lie on a line, and a point c, tells if point c is to right of oriented line AB
 */

/**
 * @brief Returns a pointer to a Grid Object on the heap containing the image's intensity matrix
 * 
 * @param filename 
 * @return Grid* 
 */
Grid* readImage(string filename){
    //open image file
    /*
    std::ifstream file("filename");
    std::vector<int> array;
    int number;
    while(file >> number) {
        array.push_back(number);
    }

    */
    ifstream imgFile(filename);
    string enc;
    imgFile >> enc;
    int width, height;
    imgFile >> width >> height;
    Grid* gr = new Grid(width,height);
    for(int r=0; r < height; r++){
        for(int c=0; c < width; c++){
            int red,blue,green;
            imgFile >> red >> blue >> green;
            int gray = (red + blue + green)/3;
            gr->setPixel(r,c,gray);
        }
    }
    imgFile.close();
    return gr;


}
/**
 * @brief applies the sobel operator to the image.
 * 
 * @param gr Pointer to grid on the heap
 * @param threshold threshold to evaluate |G| = sqrt(Gx^2 + Gy^2)
 */
void sobelOperator(Grid *gr,int threshold){
    Grid& myGrid = *gr;
    vector<vector<int>> gx = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> gy = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> xdirsobel = 
    {{-1,0,1},
    {-2,0,2},
    {-1,0,1}};
    vector<vector<int>> ydirsobel = 
    {{-1,-2,-1},
    {0,0,0},
    {1,2,1}};
    //fill gx
    for(size_t r=1; r < gx.size()-1; r++){
        for(size_t c=1; c < gx[0].size()-1; c++){
            int val = 0;
            int gxstartrow = r-1;
            int gxstartcol = c-1;
            for(size_t i=0; i < xdirsobel.size();i++){
                for(size_t j=0; j < xdirsobel.size();j++){
                    int temp =myGrid.getPixel(gxstartrow+i,gxstartcol+j)*xdirsobel[i][j];
                    val += temp;
                }
            }
            gx[r][c] = val;
        }
    }
    //fill gy
    for(size_t r=1; r < gy.size()-1; r++){
        for(size_t c=1; c < gy[0].size()-1; c++){
            int val = 0;
            int gystartrow = r-1;
            int gystartcol = c-1;
            for(size_t i=0; i < ydirsobel.size();i++){
                for(size_t j=0; j < ydirsobel.size();j++){
                    int temp =myGrid.getPixel(gystartrow+i,gystartcol+j)*ydirsobel[i][j];
                    val += temp;
                }
            }
            gy[r][c] = val;
        }
    }
    //compute magG
    vector<vector<int>> magG = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    for(size_t r=0; r < magG.size();r++){
        for(size_t c=0; c < magG[0].size();c++){
            magG[r][c] = gx[r][c]*gx[r][c] + gy[r][c] * gy[r][c];
        }
    }
    //apply threshold
    int thresholdsquared = threshold*threshold;
    Grid edges(myGrid.getWidth(),myGrid.getHeight());
    for(size_t r=0; r < magG.size();r++){
        for(size_t c=0; c < magG[0].size();c++){
            if(magG[r][c] < thresholdsquared){
                edges.setPixel(r,c,0);
            }else{
                edges.setPixel(r,c,1);
            }
        }
    }
    edges.toPPM("imagem.ppm",1);


}
void part1(){
    // readImage();
    Grid* gr = readImage("image.ppm");
    gr->toPPM("imageg.ppm",255);
    sobelOperator(gr,275);
    delete gr;
}

/**
 * @brief Performs Nonmax suprression
 * 
 * @param gr Pointer to the Grid object with the intensity matrix
 */
Grid* nonMaxSuppression(Grid* gr){
    Grid& myGrid = *gr;
    vector<vector<int>> gx = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> gy = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    vector<vector<int>> xdirsobel = 
    {{-1,0,1},
    {-2,0,2},
    {-1,0,1}};
    vector<vector<int>> ydirsobel = 
    {{-1,-2,-1},
    {0,0,0},
    {1,2,1}};
    //fill gx
    for(size_t r=1; r < gx.size()-1; r++){
        for(size_t c=1; c < gx[0].size()-1; c++){
            int val = 0;
            int gxstartrow = r-1;
            int gxstartcol = c-1;
            for(size_t i=0; i < xdirsobel.size();i++){
                for(size_t j=0; j < xdirsobel.size();j++){
                    int temp =myGrid.getPixel(gxstartrow+i,gxstartcol+j)*xdirsobel[i][j];
                    val += temp;
                }
            }
            gx[r][c] = val;
        }
    }
    //fill gy
    for(size_t r=1; r < gy.size()-1; r++){
        for(size_t c=1; c < gy[0].size()-1; c++){
            int val = 0;
            int gystartrow = r-1;
            int gystartcol = c-1;
            for(size_t i=0; i < ydirsobel.size();i++){
                for(size_t j=0; j < ydirsobel.size();j++){
                    int temp =myGrid.getPixel(gystartrow+i,gystartcol+j)*ydirsobel[i][j];
                    val += temp;
                }
            }
            gy[r][c] = val;
        }
    }
    //compute magG
    vector<vector<int>> magG = vector<vector<int>>(myGrid.getHeight(),vector<int>(myGrid.getWidth(),0));
    for(size_t r=0; r < magG.size();r++){
        for(size_t c=0; c < magG[0].size();c++){
            magG[r][c] = gx[r][c]*gx[r][c] + gy[r][c] * gy[r][c];
        }
    }
    //compute the tangent matrix
    vector<vector<double>> tangents = vector<vector<double>>(myGrid.getHeight(),vector<double>(myGrid.getWidth(),0));
    for(int r=1; r < tangents.size()-1;r++){
        for(int c=1; c < tangents[0].size()-1;c++){
            tangents[r][c] = gy[r][c]*1.0/gx[r][c];
        }
    }
    //create the suppressed grid
    Grid* suppressed = new Grid(tangents[0].size(),tangents.size());
    for(int r=1; r < suppressed->getHeight()-1;r++){
        for(int c=1; c < suppressed->getWidth()-1; c++){
            //tangent = t
            double t = tangents[r][c];
            double p1 = 0;
            double p2 = 0;
            //if t < -TAN675 (tan(67.5))
                //go vertical
            //else if t < -1
                //godiagonal
            //else if t < 0
                //go horizonal
            //else if t < 22.5:
                //go horizonal
            //else if t < 1
                // go diagonal
            //else if t < 67.5 go diagonal
            //else
                //go vertical
            if(t < -TAN675){
                p1 = myGrid.getPixel(r-1,c);
                p2 = myGrid.getPixel(r+1,c);
            }else if(t < -1){
                p1 = myGrid.getPixel(r+1,c+1);
                p2 = myGrid.getPixel(r-1,c-1);
            }else if(t < -TAN225){
                p1 = myGrid.getPixel(r+1,c+1);
                p2 = myGrid.getPixel(r-1,c-1);
            }else if(t < 0){
                p1 = myGrid.getPixel(r,c-1);
                p2 = myGrid.getPixel(r,c+1);
            }else if(t < TAN225){
                p1 = myGrid.getPixel(r,c-1);
                p2 = myGrid.getPixel(r,c+1);
            }else if(t < 1){
                p1 = myGrid.getPixel(r+1,c+1);
                p2 = myGrid.getPixel(r-1,c-1);
            }else if(t < TAN675){
                p1 = myGrid.getPixel(r+1,c+1);
                p2 = myGrid.getPixel(r-1,c-1);
            }else{
                p1 = myGrid.getPixel(r-1,c);
                p2 = myGrid.getPixel(r+1,c);
            }
            if(myGrid.getPixel(r,c)> p1 && myGrid.getPixel(r,c) > p2){
                suppressed->setPixel(r,c,myGrid.getPixel(r,c));
            }else{
                suppressed->setPixel(r,c,0);
            }
        }
    }
    suppressed->toPPM("images.ppm",255);
    return suppressed;
}
/**
 * @brief Performs Double threshold suppression and hystersis. Returns hysterisis Grid 
 * 0 - low, 1 - weak, 2 - strong
 * @param gr 
 * @param low 
 * @param high 
 */
void doubleThreshold(Grid* gr, double low, double high){
    Grid &mygrid = *gr;
    const int LOW = low*255;
    const int HIGH = high*255;
    //vector<vector<int>>* threshedPix = new vector<vector<int>>(mygrid.getHeight(),vector<int>(mygrid.getWidth(),0));
    vector<vector<int>> pix = vector<vector<int>>(mygrid.getHeight(),vector<int>(mygrid.getWidth(),0));
    for(int i=0; i < pix.size();i++){
        for(int j=0; j < pix[0].size();j++){
            int color = mygrid.getPixel(i,j);
            if(color < LOW){
                pix[i][j] = 0;
            }else if(color > HIGH){
                pix[i][j] = 2;
            }else{
                pix[i][j] = 1;
            }
        }
    }
    // 0 - low, 1 - weak, 2 - strong
    Grid* hysterisis = new Grid(gr->getWidth(),gr->getHeight());
    for(int r=1; r < hysterisis->getHeight()-1;r++){
        for(int c=1; c < hysterisis->getWidth()-1;c++){
            if(pix[r][c]==1){
                int strong = 0;
                //count the number of strong pixels within vicinith of r and c
                for(int subr = r-1; subr <= r+1; subr++){
                    for(int subc = c-1; subc <=c+1; subc++){
                        if(pix[subr][subc]==2){
                            strong += 1;
                        }
                    }
                }
                if(strong > 0){
                    hysterisis->setPixel(r,c,1);
                }else{
                    hysterisis->setPixel(r,c,0);
                } 
            }
            else if(pix[r][c]==2){
                hysterisis->setPixel(r,c,1);
            }else if(pix[r][c]==0){
                hysterisis->setPixel(r,c,0);
            }
        }
    }
    hysterisis->toPPM("imageh.ppm",1);
    delete hysterisis;


}
void part2(){
    //get gx and gy
    //calculate the angle for each pixel
    //perform non-max suppression
    Grid* gr = readImage("image.ppm");
    Grid* suprressed = nonMaxSuppression(gr);
    doubleThreshold(suprressed,0.05,0.7);
    delete gr;
    delete suprressed;
}
int main()
{
    //part1();
    part2();
}