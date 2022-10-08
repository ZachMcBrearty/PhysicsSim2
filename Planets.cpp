#include<vector>

#include "Planet.cpp"

using namespace std;

ofstream file;

double globaltime = 0;
// vector<Planet> particlesystem = {Planet(1, 1, 0.1, Vector(2,0,0)),
//             Planet(2, 1, 0.1, Vector(1,0,3)),
//             Planet(3, 1, 0.1, Vector(0,1,0)),
//             Planet(4, 1, 0.1, Vector(1,1,5))};
Planet parent = Planet(1, 100, 0.01, Vector(0,0,0));
Planet child = Planet(parent, 2, 1, 0.001, 5, 5, 0);
vector<Planet> particlesystem = {parent, child};

Vector findForce(Vector d, Planet A, Planet B) {
    return d * G * A.mass * B.mass / pow(pow(d.abs(), 2) + pow((A.radius + B.radius)/2, 2), 3/2);
}

void interact(){
    vector<float> ts(particlesystem.size(), 0);
    vector<float>::iterator t;
    int check = 0;
    vector<Planet>::iterator iter1;
    vector<Planet>::iterator iter2;
    Vector d, F;
    double localtimestep;
    do {
        check = 0;
        iter1 = particlesystem.begin();
        for (int i = 0; i < particlesystem.size(); i++) {
            iter2 = particlesystem.begin() + i + 1;
            for (int j = i+1; j < particlesystem.size(); j++) {
                d = iter1->pos - iter2->pos;
                F = findForce(d, *iter1, *iter2);
                iter1->acc -= F / iter1->mass;
                iter2->acc += F / iter2->mass;
                localtimestep = d.abs() / scale * maxtime;
                iter1->setLocal(localtimestep);
                iter2->setLocal(localtimestep);
                iter2++;
            }
            iter1++;
        }
        iter1 = particlesystem.begin();
        t = ts.begin();
        for (int i = 0; i < particlesystem.size(); i++) {
            if (*t < maxtimestep) {
                iter1->update();
                *t += iter1->timestep;
            } else {
                check += 1;
            }
            iter1++; t++;
        }
    } while (check < particlesystem.size());
}

void record() {
    double kineticEnergy = 0, gravitational = 0;
    vector<Planet>::iterator iter1 = particlesystem.begin();
    for (int i = 0; i < particlesystem.size(); i++) {
        vector<Planet>::iterator iter2 = particlesystem.begin() + i + 1;
        for (int j = i+1; j < particlesystem.size(); j++) {
            double d = (iter1->pos - iter2->pos).abs();
            gravitational -= G * iter1->mass * iter2->mass / pow(pow(d, 2) + pow((iter1->radius+iter2->radius)/2, 2), 0.5);
            iter2++;
        }
        kineticEnergy += 0.5 * iter1->mass * pow((iter1->vel).abs(), 2);
        iter1++;
    }
    file << globaltime << " # " << kineticEnergy << " , " << gravitational << " , " 
        << kineticEnergy+gravitational << " # ";
    iter1 = particlesystem.begin();
    for (int i = 0; i < particlesystem.size(); i++) {
        file << "(" << iter1->pos.x << ", " << iter1->pos.y << ", " << iter1->pos.z << ")";
        if (i != particlesystem.size()-1) file << " * ";
        iter1++;
    }
    file << endl;
}

void init() {
    file.open("PlanetsSim.txt");
}

int main() { 
    init();

    while (globaltime < maxtime) {
        record();
        interact();
        globaltime += maxtimestep;
    }
    file.close();
}