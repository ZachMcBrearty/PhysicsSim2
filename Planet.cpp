#include<iostream>

#include "Vector.cpp"

class Planet {
    public:
        Planet();
        Planet(int id, double mass_, double radius_, Vector pos_);
        Planet(Planet parent, int id, double mass_,
            double rad_, double aphelion, double semimajor, double progression);
        void update();
        void print();
        void setLocal(double localtimestep);
        int id;
        double mass, radius, timestep = maxtimestep;
        Vector pos, vel, acc;
};
Planet::Planet(){
}
Planet::Planet(int id, double mass_, double radius_, Vector pos_) {
    this->id = id;
    this->mass = mass_;
    this->radius = radius_;
    this->pos = pos_;
    this->vel = Vector(0,0,0);
    this->acc = Vector(0,0,0);
}
Planet::Planet(Planet parent, int id, double mass_,
            double rad_, double aphelion, double semimajor, double progression){
    Vector pos = VectorFromSpherical(aphelion, progression, 90);
    this->id = id;
    this->mass = mass_;
    this->radius = rad_;
    this->pos = pos;

    Vector vel = pos.normalise();
    vel = vel.rotateTheta(90);
    vel *= pow(G * parent.mass * (2/aphelion - 1/semimajor), 0.5);
    this->vel = vel + parent.vel;

    this->acc = Vector(0,0,0);
}
void Planet::update(){
    this->pos += this->vel * this->timestep;
    this->vel += this->acc * this->timestep;
    this->acc = Vector(0,0,0);
}
void Planet::print() {
    std::cout << "Planet, " << id << std::endl;
}
void Planet::setLocal(double localtimestep){
    if (localtimestep < this->timestep) {
        this->timestep = localtimestep;
    }
}