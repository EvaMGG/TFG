#ifndef _MATRIZ_H
#define _MATRIZ_H

#include "tuplas.hpp"


class Matriz4x4f // estructura que guarda una matriz de 4x4 floats
{
private:
   float coe[4][4]; // coeficientes de la matriz, por columnas

public:

   Matriz4x4f();
   float getCoe(int i, int j)const;
   void setCoe(float num, int i, int j);
   void MAT_Escalado(int x, int y, int z);
   void MAT_Trasposicion(int x, int y, int z);
   void MAT_Rotacion( const float ang_gra, const float ex, const float ey, const float ez );
   Tupla4f operator*( Tupla4f & tupla )const;
   Matriz4x4f operator*( const Matriz4x4f & der );
   Matriz4x4f operator=(const Matriz4x4f &matriz);
};

#endif
