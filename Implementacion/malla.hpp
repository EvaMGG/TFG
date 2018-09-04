#ifndef _MALLA_H
#define _MALLA_H

#include <GL/glut.h>
#include <vector>
#include "matriz.hpp"



class MallaTA
{
	unsigned num_ver;
	unsigned num_tri;
	Tupla3f* ver;
	Tupla3i* tri;
   Tupla3f* normalesCaras;
   Tupla3f* normalesVertices;

   unsigned N;

   GLint tipoDibuj;
   bool modoAjedr;
   bool caras, vertices;

public:
   MallaTA();
   MallaTA(int numv, int numt);
   void inicializarMalla(int numv, int numc);
   void extraerVerTri(std::vector<float> verticesCirculo_ply, std::vector<int> carasCirculo_ply);
   void aplicarMatriz(Matriz4x4f matriz, MallaTA malla);
   void copiarTriangulos(MallaTA malla);
   Tupla3f P2_CalcularCentro(int i);
   void dibujarObjetos();
   void cambiarCaracsDibujo(GLint tipo, bool ajedrez);
   //void normalesCaras(bool c);
   //void normalesVertices(bool v);
   void P2_CrearTriangulo(int &num_tri, int ver1, int ver2, int ver3);
   void calcularNormales();
   int P2_PreparaPerfil(Tupla3f &ver1, Tupla3f &ver2, std::vector<float> pPerfil_ply);
   void P2_AnadirExtremos(int ultimo, int num_p, Tupla3f ver1, Tupla3f ver2);
   void P2_RotarPunto(int anterior1, int ultimo, const Matriz4x4f matrizRotacion);
   void crearRotacion(std::vector<float> pPerfil_ply, int M);
   void asignarTextura(int num_p, float* distancias);

};

#endif
