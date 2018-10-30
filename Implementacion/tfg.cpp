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
#define EPSILON 0.00005

using namespace std;


class Rayo{

private :

	Tupla3f origen, direccion;

public :

	Rayo(Tupla3f o, Tupla3f d);
	Tupla3f puntoRayo(double t);
	Tupla3f getOrigen();
	Tupla3f getDireccion();
};

Rayo::Rayo(Tupla3f o, Tupla3f d){
	origen = o;
	direccion = d;
}

Tupla3f Rayo::puntoRayo(double t){
	return origen + t*direccion;
}

Tupla3f Rayo::getOrigen(){
	return origen;
}

Tupla3f Rayo::getDireccion(){
	return direccion;
}



class Superficie {

private :
	Tupla3f color;

	void intervaloInicial(Superficie &sup, Rayo rayo, double &t1, double &t2);
	double RegulaFalsi(Superficie &sup, double an, double bn, Rayo rayo);
	double NewtonRaphson(Superficie &sup, double tn, Rayo rayo);
	double derivada(Rayo rayo, double t);

public :
	
	Superficie(Tupla3f col);
	Superficie(const Superficie &sup);
	virtual double interseccion(Superficie &sup, Rayo rayo);
	virtual Tupla3f normal(Tupla3f xyz);
	virtual Tupla3f getColor();
	virtual double funcion(Tupla3f xyz) = 0;
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

double Superficie::interseccion(Superficie &sup, Rayo rayo){
	double an, bn, tn;
	int num_iteraciones = 0;
	// 1. Intervalo inicial (mediante algun método heurístico)
	intervaloInicial(sup, rayo, an, bn);
	tn = an;
	if (sup.funcion(rayo.puntoRayo(an)) == 0) return an;
	else if (sup.funcion(rayo.puntoRayo(bn)) == 0) return bn;
	else if (an == bn) return -1;
	// 2. Hasta que se cumpla el criterio de parada
	while (abs(sup.funcion(rayo.puntoRayo(tn))) > 0.1 && num_iteraciones < 100){
		// 2.1. Calcular siguiente valor de la sucesión mediante Newton-Raphson
		tn = NewtonRaphson(sup, tn, rayo);
		// 2.2. Si se sale del intervalo eludir Newton-Raphson y calcular mediante Regula-Falsi
		if (tn < an || bn < tn) tn = RegulaFalsi(sup, an, bn, rayo);
		// 2.3. Si tn es cero, hemos terminado, sino cambiar el valor de la sucesión por el extremo correspondiente del intervalo
		if (sup.funcion(rayo.puntoRayo(tn)) == 0) return tn;
		else if (sup.funcion(rayo.puntoRayo(an))*sup.funcion(rayo.puntoRayo(tn)) > 0) an = tn;
		else bn = tn;
		num_iteraciones++;
	}
	// 3. Devolver aproximación a la intersección del rayo con la superficie
	return tn;
}

void Superficie::intervaloInicial(Superficie &sup, Rayo rayo, double &t1, double &t2){
	t1 = 0;
	double incremento = 0.1;
	while ( (sup.funcion(rayo.puntoRayo(t1)) * sup.funcion(rayo.puntoRayo(t1+ incremento)) > 0) && t1 < 2) t1 = t1 + incremento;
	for (int i = 0; i < 3; i++){
		t2 = 0;
		while ( sup.funcion(rayo.puntoRayo(t2)) * sup.funcion(rayo.puntoRayo(t1)) > 0 && t2 < 2){
			t2 = t2 + incremento;
		}
		if ( sup.funcion(rayo.puntoRayo(t2)) * sup.funcion(rayo.puntoRayo(t1)) > 0 ) break;
		else incremento = incremento / 3;
	}
}

double Superficie::RegulaFalsi(Superficie &sup, double an, double bn, Rayo rayo){
	if (sup.funcion(rayo.puntoRayo(an)) == sup.funcion(rayo.puntoRayo(bn))) cout << "notANumber" << endl;
	return ( ( (an*sup.funcion(rayo.puntoRayo(bn))) - (bn*sup.funcion(rayo.puntoRayo(an))) ) / ( sup.funcion(rayo.puntoRayo(bn)) - sup.funcion(rayo.puntoRayo(an)) ) );
}

double Superficie::NewtonRaphson(Superficie &sup, double tn, Rayo rayo){
	return tn - (sup.funcion(rayo.puntoRayo(tn))/sup.derivada(rayo, tn));
}

Tupla3f Superficie::normal(Tupla3f xyz){
	return normalized( Tupla3f(funcion( xyz + Tupla3f(EPSILON, 0.0, 0.0) ) - funcion( xyz - Tupla3f(EPSILON, 0.0, 0.0) ), funcion( xyz + Tupla3f(0.0, EPSILON, 0.0) ) - funcion( xyz - Tupla3f(0.0, EPSILON, 0.0) ), funcion( xyz + Tupla3f(0.0, 0.0, EPSILON) ) - funcion( xyz - Tupla3f(0.0, 0.0, EPSILON) )) );
}

double Superficie::derivada(Rayo rayo, double t){
	return (funcion(rayo.puntoRayo(t + EPSILON)) + funcion(rayo.puntoRayo(t))) / EPSILON;
}








class Elipse : public Superficie{

private:
	double radio0, radio1, radio2;

public :

	Elipse(double rad0, double rad1, double rad2, Tupla3f color);
	Elipse(const Elipse &eli);
	double funcion(Tupla3f xyz);
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

double Elipse::funcion(Tupla3f xyz){
	return ((pow(xyz.coo[0],2) / pow(radio0,2)) + (pow(xyz.coo[1],2) / pow(radio1,2)) + (pow(xyz.coo[2],2) / pow(radio2,2)) -1);
}





class OvaloidePrueba : public Superficie{

private:


public :

	OvaloidePrueba(Tupla3f c);
	double funcion(Tupla3f xyz);
};



OvaloidePrueba::OvaloidePrueba(Tupla3f c):Superficie(c){

}

double OvaloidePrueba::funcion(Tupla3f xyz){
	return (pow(xyz.coo[0],4) + pow(xyz.coo[1],4) + pow(xyz.coo[2],4) -3);
}



class Cuadrica : public Superficie{

private:

	double a, b, c, d, e, f, g, h, i, j;

public :

	Cuadrica(double ac, double bc, double cc, double dc, double ec, double fc, double gc, double hc, double ic, double jc, Tupla3f c);
	double funcion(Tupla3f xyz);
};



Cuadrica::Cuadrica(double ac, double bc, double cc, double dc, double ec, double fc, double gc, double hc, double ic, double jc, Tupla3f col):Superficie(col){
	a = ac;
	b = bc;
	c = cc;
	d = dc;
	e = ec;
	f = fc;
	g = gc;
	h = hc;
	i = ic;
	j = jc;
}

double Cuadrica::funcion(Tupla3f xyz){
	return (a*pow(xyz.coo[0],2) + b*(xyz.coo[0]*xyz.coo[1]) + c*(xyz.coo[0]*xyz.coo[2]) +d*xyz.coo[0] +e*pow(xyz.coo[1],2) +f*(xyz.coo[1]*xyz.coo[2]) +g*xyz.coo[2] + h*pow(xyz.coo[2],2) + i*xyz.coo[2] +j);
}







class Iluminacion{
private:
	Tupla3f color;

public:
	Iluminacion(Tupla3f col);
	virtual Tupla3f LeyLambert(Tupla3f objColor, Tupla3f normal, Tupla3f punto);
	virtual Tupla3f getDireccion(Tupla3f punto)=0;
};

Iluminacion::Iluminacion(Tupla3f col){
	color = col;
}

Tupla3f Iluminacion::LeyLambert(Tupla3f objColor, Tupla3f normal, Tupla3f punto) {
	Tupla3f L = Tupla3f(color.coo[0]*objColor.coo[0], color.coo[1]*objColor.coo[1], color.coo[2]*objColor.coo[2]) * (normal | getDireccion(punto));
	return L;
}



class LuzDireccional : public Iluminacion {
private:
	Tupla3f direccion;

public:
	LuzDireccional(Tupla3f d, Tupla3f c);
	Tupla3f getDireccion(Tupla3f punto);
};

LuzDireccional::LuzDireccional(Tupla3f d, Tupla3f c):Iluminacion(c) {
	direccion = d;
}

Tupla3f LuzDireccional::getDireccion(Tupla3f punto) {
	return direccion;
}


class LuzPuntual : public Iluminacion {
private:
	Tupla3f posicion;
	Tupla3f color;

public:
	LuzPuntual(Tupla3f p, Tupla3f c);
	Tupla3f getDireccion(Tupla3f punto);
};

LuzPuntual::LuzPuntual(Tupla3f p, Tupla3f c):Iluminacion(c) {
	posicion = p;
	color = c;
}

Tupla3f LuzPuntual::getDireccion(Tupla3f punto) {
	return normalized(posicion - punto);
}






vector<Superficie *> superficies;
vector<Iluminacion *> fuentesLuz;

double a1, a2, b1, b2;
double s = 6.0;
Tupla3f e = Tupla3f(0,0,8);
Tupla3f g = Tupla3f(0,0,-2);
Tupla3f vup = Tupla3f(0,1,0);
Tupla3f u, v, w;
Tupla3f *image;



void Inicializar(int x, int y) {

	//Inicialización de los superficies de la escena

	/*superficies.push_back(new Esfera(Tupla3f(0,0,0), 2.0, Tupla3f(0.3, 0.3, 0.3)));
	superficies.push_back(new Cubo(Tupla3f(1.5,1.5,1.5), Tupla3f(2,2,2), Tupla3f(0.5, 0.5, 0.5)));
	superficies.push_back(new Esfera(Tupla3f(-1.5,-1.5,1.5), 0.5, Tupla3f(0.55, 0.55, 0.55)));*/

	//superficies.push_back(new Elipse(2, 1.6, 1.3, Tupla3f(0.5, 0.1, 0.1)));
	superficies.push_back(new OvaloidePrueba(Tupla3f(0.1, 0.3, 0.2)));

	//cono
	//superficies.push_back(new Cuadrica(0.5,0.0,0.0,0.0,-0.5,0.0,0.0,0.5,0.0,0.0, Tupla3f(0.5, 0.1, 0.1)));

	//hiperboloide de una hoja
	//superficies.push_back(new Cuadrica(0.5,0.0,0.0,0.0,-0.5,0.0,0.0,-0.5,0.0,-0.5, Tupla3f(0.5, 0.1, 0.1)));

	fuentesLuz.push_back(new LuzDireccional(Tupla3f(1.5,1.5,1.65), Tupla3f(1, 1, 1)));
	fuentesLuz.push_back(new LuzPuntual(Tupla3f(-2,-2,2), Tupla3f(1, 1, 1)));




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
void primeraInterseccion(Rayo rayo, double &min_inter, int &indice_min){
	indice_min = 0;
	min_inter = (superficies[0])->interseccion(*(superficies[0]), rayo);

	for (int k = 1; k < (int)superficies.size(); k++) {
		double intersec = (superficies[k])->interseccion(*(superficies[k]), rayo);

		if ( (intersec < min_inter && intersec >= 0) || (min_inter == -1 && intersec >= 0) ) {
			indice_min = k;
			min_inter = intersec;
		}
	}
}

//Devuelve true si hay una superficie que se interpone entre la superficie cuyo indice se pasa como parámetro y la luz
bool interposicionSuperficie(Rayo rayo, double min_inter, int indice_min){
	bool inter = false;
	int k=0;
	while ( (k < (int)superficies.size()) && !inter) {
		if ( ((superficies[k])->interseccion(*(superficies[k]),  Rayo(rayo.puntoRayo(min_inter), fuentesLuz[0]->getDireccion(rayo.puntoRayo(min_inter))) ) >= 0.0) && (k!= indice_min) )  inter = true;
		k++;
	}
	return inter;
}


//ARREGLAR CONSTRUCTOR DE COPIA









Tupla3f* Image(int x, int y){
	//Dirección de los rayos de visión y matriz que contiene a la imagen
	Tupla3f direccion;
	double int_ant;

	for(int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++){
			direccion = (a1 + (b1 - a1)*(i/(x-1.0)))*u + (a2 + (b2 - a2)*(j/(y-1.0)))*v -s*w;
			Rayo rayo = Rayo(e, direccion);

			if (!superficies.empty()){

				int indice_min;

				primeraInterseccion(rayo, int_ant, indice_min);

				if (int_ant >= 0) {
					image[i*y +j] = fuentesLuz[0]->LeyLambert( (superficies[indice_min])->getColor(), (superficies[indice_min])->normal(rayo.puntoRayo(int_ant)), rayo.puntoRayo(int_ant));

					if (interposicionSuperficie(rayo, int_ant, indice_min) == true) image[i*y +j] = Tupla3f(0, 0, 0);
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









