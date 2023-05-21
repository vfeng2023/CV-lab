#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;
int main(){
    //cout << "Hello world!";
    ofstream myFile;
    myFile.open("test.txt");
    for(int i=0; i < 10; i++)
        myFile << setprecision(15) << 4.565675467546745675467456746846756747;
    myFile.close();

}