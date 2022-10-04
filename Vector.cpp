#include<math.h>
#include "Constants.cpp"
class Vector {
    public:
        Vector();
        Vector(double x, double y);
        Vector(double x, double y, double z);
        double x, y, z;
        Vector rotateTheta(double theta);
        Vector operator + (Vector v);
        void operator += (Vector v);
        Vector operator - (Vector v);
        void operator -= (Vector v);
        Vector operator * (double q);
        void operator *= (double q);
        Vector operator / (double q);
        double abs();
        Vector normalise();
        void print();
        void print(char* end);
        
};
Vector::Vector(){

}
Vector::Vector(double x_, double y_) {
    this->x = x_; this->y=y_; this->z=0;
}

Vector::Vector(double x_, double y_, double z_) {
    this->x = x_; this->y=y_; this->z=z_;
}

Vector Vector::rotateTheta(double theta) {
    // theta radians
    double x_ = this->x*cos(theta) - this->y*sin(theta);
    double y_ = this->x*sin(theta) + this->y*cos(theta);
    return Vector(x_, y_, this->z);
}

Vector Vector::operator +(Vector v) {
    return Vector(this->x+v.x, this->y+v.y, this->z+v.z);
}
void Vector::operator +=(Vector v) {
    this->x+=v.x;
    this->y+=v.y;
    this->z+=v.z;
}
Vector Vector::operator -(Vector v) {
    return Vector(this->x-v.x, this->y-v.y, this->z-v.z);
}
void Vector::operator -=(Vector v) {
    this->x-=v.x;
    this->y-=v.y;
    this->z-=v.z;
}
Vector Vector::operator *(double q) {
    return Vector(this->x*q, this->y*q, this->z*q);
}
void Vector::operator *= (double q){
    this->x *= q;
    this->y *= q;
    this->z *= q;
}
Vector Vector::operator /(double q) {
    return Vector(this->x/q, this->y/q, this->z/q);
}
double Vector::abs() {
    return pow(pow(this->x,2)+pow(this->y, 2)+pow(this->z,2), 0.5);
}
Vector Vector::normalise() {
    return *this / this->abs();
}
void Vector::print(){
    std::cout << "(" << this->x << ", " << this->y << ", " << this->z << ")" << std::endl;
}
void Vector::print(char* end){
    std::cout << "(" << this->x << ", " << this->y << ", " << this->z << ")" << end;
}

Vector VectorFromSpherical(double r, double theta, double phi) {
    double x = r*cos(theta)*sin(phi);
    double y = r*sin(theta)*sin(phi);
    double z = r*cos(phi);
    return Vector(x, y, z);
}