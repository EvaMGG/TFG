// *********************************************************************
// **
// ** Informática Gráfica, curso 2014-15
// ** 
// ** Práctica 1  (implementación)
// **
// *********************************************************************


#include <GL/glut.h>
#include <iostream>
#include "error-ogl.hpp"
#include "tuplas.hpp"   // Tupla3f, Tupla4f, Tupla3i
#include "tfg.hpp"
//#include "objeto.hpp"
//#include "esfera.hpp"
//#include "cubo.hpp"


#include <vector>
#include "file_ply_stl.hpp"
#include <string>
#include <fstream>  // ifstream
#include <vector>


#define INFINITO 99999 

using namespace std;



class Objeto {

private :


public :
	
	Objeto();
	Objeto(const Objeto &obj);
	virtual double interseccion(Tupla3f origen, Tupla3f direccion) = 0;
	virtual Tupla3f getColor()=0;
	virtual Tupla3f normal(Tupla3f e, Tupla3f d, float t)=0;

};

Objeto::Objeto(){}
Objeto::Objeto(const Objeto &obj){}



class Esfera : public Objeto {

private :

	Tupla3f centro;
	float radio;
	Tupla3f color;

public :

	Esfera(Tupla3f c, float r, Tupla3f col);
	Esfera(const Esfera &sph);
	double interseccion(Tupla3f origen, Tupla3f direccion);
	Tupla3f getColor();
	Tupla3f normal(Tupla3f e, Tupla3f d, float t);

};

Esfera::Esfera(Tupla3f c, float r, Tupla3f col):Objeto() {
	centro = c;
	radio = r;
	color = col;
}


Esfera::Esfera(const Esfera &sph): Objeto(sph) {
	centro = sph.centro;
	radio = sph.radio;
	color = sph.color;
}


double Esfera::interseccion(Tupla3f o, Tupla3f d){
	double interseccion = -1;

	float discriminante = (d|(o-centro))*(d|(o-centro)) - (d|d)*(((o-centro)|(o-centro)) - radio*radio);

	if (discriminante > 0) {
		

		double t0 = ( ( -(d|(o-centro)) + sqrt( discriminante ) ) / ( d|d ) );
		double t1 = ( ( -(d|(o-centro)) - sqrt( discriminante ) ) / ( d|d ) );

		if (t0 > 0) interseccion = t0;
		if (t0 > t1 && t1>0) interseccion = t1;
	}
	return interseccion;
}

Tupla3f Esfera::getColor() {
	return color;
}

Tupla3f Esfera::normal(Tupla3f e, Tupla3f d, float t) {
	return Tupla3f( ((e + t*d)-centro)/radio );
}



class Cubo : public Objeto {

private :

	Tupla3f esquina1, esquina2, color;

public :

	Cubo(Tupla3f e1, Tupla3f e2, Tupla3f e3);
	Cubo(const Cubo &squ);
	double interseccion(Tupla3f origen, Tupla3f direccion);
	Tupla3f getColor();
	Tupla3f normal(Tupla3f e, Tupla3f d, float t);

};




Cubo::Cubo(Tupla3f e1, Tupla3f e2, Tupla3f e3):Objeto() {
	esquina1 = e1;
	esquina2 = e2;
	color = e3;
}

Cubo::Cubo(const Cubo &squ):Objeto(squ) {
	esquina1 = squ.esquina1;
	esquina2 = squ.esquina2;
	color = squ.color;
}

double Cubo::interseccion(Tupla3f origen, Tupla3f direccion){
	float tx_min, tx_max, ty_min, ty_max, tz_min, tz_max, t0, t1;


	if (direccion.coo[0] > 0) {
		tx_min = (esquina1.coo[0] - origen.coo[0])/direccion.coo[0];
		tx_max = (esquina2.coo[0] - origen.coo[0])/direccion.coo[0];
	}
	else {
		tx_min = (esquina2.coo[0] - origen.coo[0])/direccion.coo[0];
		tx_max = (esquina1.coo[0] - origen.coo[0])/direccion.coo[0];
	}

	if (direccion.coo[1] > 0) {
		ty_min = (esquina1.coo[1] - origen.coo[1])/direccion.coo[1];
		ty_max = (esquina2.coo[1] - origen.coo[1])/direccion.coo[1];
	}
	else {
		ty_min = (esquina2.coo[1] - origen.coo[1])/direccion.coo[1];
		ty_max = (esquina1.coo[1] - origen.coo[1])/direccion.coo[1];
	}

	if (direccion.coo[2] > 0) {
		tz_min = (esquina1.coo[2] - origen.coo[2])/direccion.coo[2];
		tz_max = (esquina2.coo[2] - origen.coo[2])/direccion.coo[2];
	}
	else {
		tz_min = (esquina2.coo[2] - origen.coo[2])/direccion.coo[2];
		tz_max = (esquina1.coo[2] - origen.coo[2])/direccion.coo[2];
	}

	if (tx_min > ty_min) t0 = tx_min;
	else t0 = ty_min;
	if (tz_min > t0) t0 = tz_min;

	if (tx_max < ty_max) t1 = tx_max;
	else t1 = ty_max;
	if (tz_max < t1) t1 = tz_max;


	if (t0 < t1) return t0;
	else return -1;
}

Tupla3f Cubo::getColor() {
	return color;
}

Tupla3f Cubo::normal(Tupla3f e, Tupla3f d, float t) {
	if ( abs((e+t*d).coo[0] - esquina2.coo[0]) <= 0.01 ) return Tupla3f(1,0,0);
	else if ( abs((e+t*d).coo[1] - esquina2.coo[1]) <= 0.01 ) return Tupla3f(0,1,0);
	else if ( abs((e+t*d).coo[2] - esquina2.coo[2]) <= 0.01 ) return Tupla3f(0,0,1);
	else if ( abs((e+t*d).coo[0] - esquina1.coo[0]) <= 0.01 ) return Tupla3f(-1,0,0);
	else if ( abs((e+t*d).coo[1] - esquina1.coo[1]) <= 0.01 ) return Tupla3f(0,-1,0);
	else if ( abs((e+t*d).coo[2] - esquina1.coo[2]) <= 0.01 ) return Tupla3f(0,0,-1);
}







class LuzDireccional{
private:
	Tupla3f direccion;
	Tupla3f color;

public:
	LuzDireccional(Tupla3f d, Tupla3f c);
	Tupla3f LeyLambert(Tupla3f objColor, Tupla3f normal);
	Tupla3f getDireccion();
};

LuzDireccional::LuzDireccional(Tupla3f d, Tupla3f c) {
	direccion = d;
	color = c;
}

Tupla3f LuzDireccional::LeyLambert(Tupla3f objColor, Tupla3f normal) {
	Tupla3f L = Tupla3f(color.coo[0]*objColor.coo[0], color.coo[1]*objColor.coo[1], color.coo[2]*objColor.coo[2]) * (normal | direccion);
	return L;
}

Tupla3f LuzDireccional::getDireccion() {
	return direccion;
}


class LuzPuntual {
private:
	Tupla3f posicion;
	Tupla3f color;

public:
	LuzPuntual(Tupla3f p, Tupla3f c);
	Tupla3f LeyLambert(Tupla3f objColor, Tupla3f normal, Tupla3f punto);
	Tupla3f getDireccion(Tupla3f punto);
};

LuzPuntual::LuzPuntual(Tupla3f p, Tupla3f c) {
	posicion = p;
	color = c;
}

Tupla3f LuzPuntual::LeyLambert(Tupla3f objColor, Tupla3f normal, Tupla3f punto) {
	Tupla3f L = Tupla3f(color.coo[0]*objColor.coo[0], color.coo[1]*objColor.coo[1], color.coo[2]*objColor.coo[2]) * (normal | LuzPuntual::getDireccion(punto));
	return L;
}

Tupla3f LuzPuntual::getDireccion(Tupla3f punto) {
	return normalized(posicion - punto);
}



vector<Objeto *> objetos;

float a1, a2, b1, b2;
float s = 6.0;
Tupla3f e = Tupla3f(0,0,8);
Tupla3f g = Tupla3f(0,0,-2);
Tupla3f vup = Tupla3f(0,1,0);
Tupla3f u, v, w;
Tupla3f *image;

LuzDireccional luz = LuzDireccional(Tupla3f(1.5,1.5,1.65), Tupla3f(1, 1, 1));


void TFG_Inicializar(int x, int y) {

	//Inicialización de los objetos de la escena

	objetos.push_back(new Esfera(Tupla3f(0,0,0), 2.0, Tupla3f(0.3, 0.3, 0.3)));
	objetos.push_back(new Cubo(Tupla3f(1.5,1.5,1.5), Tupla3f(2,2,2), Tupla3f(0.5, 0.5, 0.5)));
	objetos.push_back(new Esfera(Tupla3f(-1.5,-1.5,1.5), 0.5, Tupla3f(0.55, 0.55, 0.55)));

	image = new Tupla3f [x*y];

	//Inicialización de los elementos de la vista

	//Nos aseguramos de que la imagen no aparezca ensanchada
	a1=-x/200.0;
	a2=y/200.0;
	b1=-a1;
	b2=-a2;

	//Calculamos la base de visualización
	w = (-1)*normalized(g);
	u = normalized(vup*w);
	v = w*u;
}

//Devuelve en los parametros pasados por referencia el valor de t en el que se corta por primera vez un objeto y el indice que indica la posición en el vector de objetos del objeto al que el rayo interseca por primera vez.
void primeraInterseccion(Tupla3f direccion, double &min_inter, int &indice_min){
	indice_min = 0;
	min_inter = (objetos[0])->interseccion(e,direccion);

	for (int k = 1; k < (int)objetos.size(); k++) {
		double intersec = (objetos[k])->interseccion(e,direccion);

		if ( (intersec < min_inter && intersec >= 0) || (min_inter == -1 && intersec >= 0) ) {
			indice_min = k;
			min_inter = intersec;
		}
	}
}

//Devuelve true si hay un objeto que se interpone entre el objeto cuyo indice se pasa como parámetro y la luz
bool interposicionObjeto(Tupla3f direccion, double min_inter, int indice_min){
	bool inter = false;
	int k=0;
	while ( (k < (int)objetos.size()) && !inter) {
		if ( ((objetos[k])->interseccion( e + min_inter*direccion, luz.getDireccion() ) >= 0.0) && (k!= indice_min) )  inter = true;
		k++;
	}
	return inter;
}

Tupla3f* TFG_Image(int x, int y){
	//Dirección de los rayos de visión y matriz que contiene a la imagen
	Tupla3f direccion;
	//Tupla3f *image = new Tupla3f [x*y];

	for(int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++){
			direccion = (a1 + (b1 - a1)*(i/(x-1.0)))*u + (a2 + (b2 - a2)*(j/(y-1.0)))*v -s*w;
			

			if (!objetos.empty()){

				int indice_min;
				double int_ant;

				primeraInterseccion(direccion, int_ant, indice_min);

				if (int_ant >= 0) {
					image[i*y +j] = luz.LeyLambert( (objetos[indice_min])->getColor(), (objetos[indice_min])->normal(e, direccion, int_ant));
					//std::cout << ((objetos[indice_min])->normal(e, direccion, int_ant))).coo[0] << endl;
					if (interposicionObjeto(direccion, int_ant, indice_min) == true) image[i*y +j] = Tupla3f(0, 0, 0);
				}
				else image[i*y +j] = Tupla3f(1, 1, 1);
			}
		}
    	}
	return image;
}


void TFG_Destruir()
{
	for (int i = 0; i < (int)objetos.size(); i++){
		delete objetos[i];
	}
	delete image;
}









