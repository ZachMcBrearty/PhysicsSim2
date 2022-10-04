#include<iostream>
#include<math.h>
#include<vector>
#include "Planet.cpp"

using namespace std;

Vector findForce(Planet A, Planet B) {
    Vector d = A.pos - B.pos;
    return d * G * A.mass * B.mass / pow(pow(d.abs(), 2) + pow((A.radius + B.radius)/2, 2), 3/2);
}

void interact(vector<Planet> system){
    vector<Planet>::iterator iter1 = system.begin();
    for (int i = 0; i < system.size(); i++) {
        vector<Planet>::iterator iter2 = system.begin() + i + 1;
        for (int j = i+1; j < system.size(); j++) {
            Vector F = findForce(*iter1, *iter2);
            iter1->acc -= F / iter1->mass;
            iter2->acc += F / iter2->mass;
            iter1->print();
            iter2->print();
            iter2++;
        }
        iter1++;
    } 
}

void init() {

}

int main() {
    vector<Planet> system = {Planet(1, 0, 0, Vector(0,0,0)),
                             Planet(2, 0, 0, Vector(0,0,0)),
                             Planet(3, 0, 0, Vector(0,0,0)),
                             Planet(4, 0, 0, Vector(0,0,0))};
    interact(system);
} 