#include<math.h>
#include<iostream>
#include<fstream>

const double G = 1; // 6.67 * pow(10, -11);
const long long scale = 10; //5 * pow(10, 18); //2 * pow(10, 11)
const int maxtimestep = 0.001;
const int maxtime = 10;
int stepcount = 10000 + 1;