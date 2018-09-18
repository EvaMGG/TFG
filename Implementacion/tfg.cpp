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
//#include "superficie.hpp"
//#include "esfera.hpp"
//#include "cubo.hpp"


#include <vector>
#include "file_ply_stl.hpp"
#include <string>
#include <fstream>  // ifstream
#include <vector>
#include <math.h>


#define INFINITO 99999 

using namespace std;



class Superficie {

private :
	Tupla3f color;

public :
	
	Superficie(Tupla3f col);
	Superficie(const Superficie &sup);
	virtual double interseccion(Tupla3f origen, Tupla3f direccion) = 0;
	virtual Tupla3f getColor();
	virtual Tupla3f normal(Tupla3f e, Tupla3f d, double t) = 0;
	virtual double funcion(Tupla3f o, Tupla3f d, double t) = 0;
	virtual double derivada(Tupla3f o, Tupla3f d, double t) = 0;

};

Superficie::Superficie(Tupla3f col){
	color = col;
}

Superficie::Superficie(const Superficie &sup){
	color = sup.color;
}


Tupla3f Superficie::getColor(){
	return color;
}








class Elipse : public Superficie{

private:
	double radio0, radio1, radio2;

	void intervaloInicial(Superficie &sup, Tupla3f o, Tupla3f d, double &t1, double &t2);
	double RegulaFalsi(Superficie &sup, double an, double bn, Tupla3f o, Tupla3f d);
	double NewtonRaphson(Superficie &sup, double tn, Tupla3f o, Tupla3f d);

public :

	Elipse(double rad0, double rad1, double rad2, Tupla3f color);
	Elipse(const Elipse &eli);
	double interseccion(Superficie &sup, Tupla3f o, Tupla3f d);
	Tupla3f normal(Tupla3f e, Tupla3f d, double t);
	double funcion(Tupla3f o, Tupla3f d, double t);
	double derivada(Tupla3f o, Tupla3f d, double t);
};



Elipse::Elipse(double rad0, double rad1, double rad2, Tupla3f col):Superficie(col){
	radio0 = rad0;
	radio1 = rad1;
	radio2 = rad2;
}

Elipse::Elipse(const Elipse &eli):Superficie(eli){
	radio0 = eli.radio0;
	radio1 = eli.radio1;
	radio2 = eli.radio2;
}



void Elipse::intervaloInicial(Superficie &sup, Tupla3f o, Tupla3f d, double &t1, double &t2){
	t1 = 0;
	double incremento = 0.25;
	for (int i = 0; i < 5; i++){
		t2 = 0;
		while ( sup.funcion(o, d, t2) * sup.funcion(o, d, t1) > 0 && t2 < 7){
			t2 = t2 + incremento;
		}
		if ( sup.funcion(o, d, t2) * sup.funcion(o, d, t1) > 0 ) break;
		else incremento = incremento / 3;
	}
}

double Elipse::RegulaFalsi(Superficie &sup, double an, double bn, Tupla3f o, Tupla3f d){
	return ( ( (an*sup.funcion(o,d,bn)) - (bn*sup.funcion(o,d,an)) ) / ( sup.funcion(o,d,bn) - sup.funcion(o,d,an) ) );
}

double Elipse::NewtonRaphson(Superficie &sup, double tn, Tupla3f o, Tupla3f d){
	return tn - (sup.funcion(o, d, tn)/sup.derivada(o, d, tn));
}

/*double Elipse::AlgoritmoInterseccion(Superficie &sup, Tupla3f o, Tupla3f d){

}*/




double Elipse::interseccion(Superficie &sup, Tupla3f o, Tupla3f d){
	double an, bn, tn=0;
	// 1. Intervalo inicial (mediante algun método heurístico)
	intervaloInicial(sup, o, d, an,bn);
	// 2. Hasta que se cumpla el criterio de parada
	while (abs(sup.funcion(o, d, tn)) > 0.01){
		// 2.1. Calcular siguiente valor de la sucesión mediante Newton-Raphson
		tn = NewtonRaphson(sup, tn, o, d);
		// 2.2. Si se sale del intervalo eludir Newton-Raphson y calcular mediante Regula-Falsi
		if (tn < an || bn < tn) tn = RegulaFalsi(sup, an, bn, o, d);
		// 2.3. Si tn es cero, hemos terminado, sino cambiar el valor de la sucesión por el extremo correspondiente del intervalo
		if (sup.funcion(o, d, tn) == 0) return tn;
		else if (an*tn > 0) an = tn;
		else bn = tn;
	}
	// 3. Devolver aproximación a la intersección del rayo con la superficie
	return tn;
}

Tupla3f Elipse::normal(Tupla3f e, Tupla3f d, double t){
	return Tupla3f(1.0,1.0,1.0);
}




double Elipse::funcion(Tupla3f o, Tupla3f d, double t){
	return ((pow(o.coo[0]+ d.coo[0]*t,2) / pow(radio0,2)) + (pow(o.coo[1]+ d.coo[1]*t,2) / pow(radio1,2)) + (pow(o.coo[2]+ d.coo[2]*t,2) / pow(radio2,2)) -1);
}

double Elipse::derivada(Tupla3f o, Tupla3f d, double t){
	return 2*((o.coo[0]+ d.coo[0]*t / pow(radio0,2)) + (o.coo[1]+ d.coo[1]*t / pow(radio1,2)) + (o.coo[2]+ d.coo[2]*t / pow(radio2,2)) );
}
















class Esfera : public Superficie {

private :

	Tupla3f centro;
	double radio;

public :

	Esfera(Tupla3f c, double r, Tupla3f col);
	Esfera(const Esfera &sph);
	double interseccion(Tupla3f origen, Tupla3f direccion);
	Tupla3f normal(Tupla3f e, Tupla3f d, double t);

// BORRAR MAS TARDE
	double funcion(Tupla3f o, Tupla3f d, double t);
	double derivada(Tupla3f o, Tupla3f d, double t);

};

//BORRAR MAS TARDE
double Esfera::funcion(Tupla3f o, Tupla3f d, double t){
}
double Esfera::derivada(Tupla3f o, Tupla3f d, double t){
}







Esfera::Esfera(Tupla3f c, double r, Tupla3f col):Superficie(col) {
	centro = c;
	radio = r;
}


Esfera::Esfera(const Esfera &sph): Superficie(sph) {
	centro = sph.centro;
	radio = sph.radio;
}


double Esfera::interseccion(Tupla3f o, Tupla3f d){
	double interseccion = -1;

	double discriminante = (d|(o-centro))*(d|(o-centro)) - (d|d)*(((o-centro)|(o-centro)) - radio*radio);

	if (discriminante > 0) {
		

		double t0 = ( ( -(d|(o-centro)) + sqrt( discriminante ) ) / ( d|d ) );
		double t1 = ( ( -(d|(o-centro)) - sqrt( discriminante ) ) / ( d|d ) );

		if (t0 > 0) interseccion = t0;
		if (t0 > t1 && t1>0) interseccion = t1;
	}
	return interseccion;
}


Tupla3f Esfera::normal(Tupla3f e, Tupla3f d, double t) {
	return Tupla3f( ((e + t*d)-centro)/radio );
}



class Cubo : public Superficie {

private :

	Tupla3f esquina1, esquina2;

public :

	Cubo(Tupla3f e1, Tupla3f e2, Tupla3f e3);
	Cubo(const Cubo &squ);
	double interseccion(Tupla3f origen, Tupla3f direccion);
	Tupla3f normal(Tupla3f e, Tupla3f d, double t);

// BORRAR MAS TARDE
	double funcion(Tupla3f o, Tupla3f d, double t);
	double derivada(Tupla3f o, Tupla3f d, double t);

};

//BORRAR MAS TARDE
double Cubo::funcion(Tupla3f o, Tupla3f d, double t){
}
double Cubo::derivada(Tupla3f o, Tupla3f d, double t){
}







Cubo::Cubo(Tupla3f e1, Tupla3f e2, Tupla3f e3):Superficie(e3) {
	esquina1 = e1;
	esquina2 = e2;
}

Cubo::Cubo(const Cubo &squ):Superficie(squ) {
	esquina1 = squ.esquina1;
	esquina2 = squ.esquina2;
}

double Cubo::interseccion(Tupla3f origen, Tupla3f direccion){
	double tx_min, tx_max, ty_min, ty_max, tz_min, tz_max, t0, t1;


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


Tupla3f Cubo::normal(Tupla3f e, Tupla3f d, double t) {
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



vector<Superficie *> superficies;

double a1, a2, b1, b2;
double s = 6.0;
Tupla3f e = Tupla3f(0,0,8);
Tupla3f g = Tupla3f(0,0,-2);
Tupla3f vup = Tupla3f(0,1,0);
Tupla3f u, v, w;
Tupla3f *image;

LuzDireccional luz = LuzDireccional(Tupla3f(1.5,1.5,1.65), Tupla3f(1, 1, 1));


void Inicializar(int x, int y) {

	//Inicialización de los superficies de la escena

	/*superficies.push_back(new Esfera(Tupla3f(0,0,0), 2.0, Tupla3f(0.3, 0.3, 0.3)));
	superficies.push_back(new Cubo(Tupla3f(1.5,1.5,1.5), Tupla3f(2,2,2), Tupla3f(0.5, 0.5, 0.5)));
	superficies.push_back(new Esfera(Tupla3f(-1.5,-1.5,1.5), 0.5, Tupla3f(0.55, 0.55, 0.55)));*/

	superficies.push_back(new Elipse(0.1, 0.3, 0.5, Tupla3f(0.3, 0.3, 0.3)));




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

//Devuelve en los parametros pasados por referencia el valor de t en el que se corta por primera vez una superficie y el indice que indica la posición en el vector de superficies de la superficie a la que el rayo interseca por primera vez.
void primeraInterseccion(Tupla3f direccion, double &min_inter, int &indice_min){
	indice_min = 0;
	min_inter = (superficies[0])->interseccion(e,direccion);

	for (int k = 1; k < (int)superficies.size(); k++) {
		double intersec = (superficies[k])->interseccion(e,direccion);

		if ( (intersec < min_inter && intersec >= 0) || (min_inter == -1 && intersec >= 0) ) {
			indice_min = k;
			min_inter = intersec;
		}
	}
}

//Devuelve true si hay una superficie que se interpone entre la superficie cuyo indice se pasa como parámetro y la luz
bool interposicionSuperficie(Tupla3f direccion, double min_inter, int indice_min){
	bool inter = false;
	int k=0;
	while ( (k < (int)superficies.size()) && !inter) {
		if ( ((superficies[k])->interseccion( e + min_inter*direccion, luz.getDireccion() ) >= 0.0) && (k!= indice_min) )  inter = true;
		k++;
	}
	return inter;
}


//ARREGLAR CONSTRUCTOR DE COPIA









Tupla3f* Image(int x, int y){
	//Dirección de los rayos de visión y matriz que contiene a la imagen
	Tupla3f direccion;

	for(int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++){
			direccion = (a1 + (b1 - a1)*(i/(x-1.0)))*u + (a2 + (b2 - a2)*(j/(y-1.0)))*v -s*w;
			

			if (!superficies.empty()){

				int indice_min;
				double int_ant;

				primeraInterseccion(direccion, int_ant, indice_min);

				if (int_ant >= 0) {
					image[i*y +j] = luz.LeyLambert( (superficies[indice_min])->getColor(), (superficies[indice_min])->normal(e, direccion, int_ant));

					if (interposicionSuperficie(direccion, int_ant, indice_min) == true) image[i*y +j] = Tupla3f(0, 0, 0);
				}
				else image[i*y +j] = Tupla3f(1, 1, 1);
			}
		}
    	}
	return image;
}


void DestruirTFG()
{
	for (int i = 0; i < (int)superficies.size(); i++){
		delete superficies[i];
	}
	delete image;
}









