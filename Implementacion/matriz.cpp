#include <cmath>
#include "matriz.hpp"



Matriz4x4f::Matriz4x4f(){
   for (int i=0; i<4; i++){
      for (int j=0; j<4; j++){
         if (i==j) coe[i][j]=1;
         else coe[i][j]=0;
      }
   }
}

float Matriz4x4f::getCoe(int i, int j)const{
   return coe[i][j];
}

void Matriz4x4f::setCoe(float num, int i, int j){
   coe[i][j]=num;
}

void Matriz4x4f::MAT_Escalado(int x, int y, int z)
{
   coe[0][0]=x;
   coe[1][1]=y;
   coe[2][2]=z;
}

void Matriz4x4f::MAT_Trasposicion(int x, int y, int z)
{
   coe[0][3]=x;
   coe[1][3]=y;
   coe[2][3]=z;
}

void Matriz4x4f::MAT_Rotacion( const float ang_gra, const float ex, const float ey, const float ez ) {

   coe[0][0]=((1-cos(ang_gra))*ex)*ex+cos(ang_gra);
   coe[0][1]=((1-cos(ang_gra))*ex)*ey-sin(ang_gra)*ez;
   coe[0][2]=((1-cos(ang_gra))*ex)*ez+sin(ang_gra)*ey;

   coe[1][0]=((1-cos(ang_gra))*ey)*ex+sin(ang_gra)*ez;
   coe[1][1]=((1-cos(ang_gra))*ey)*ey+cos(ang_gra);
   coe[1][2]=((1-cos(ang_gra))*ey)*ez-sin(ang_gra)*ex;

   coe[2][0]=((1-cos(ang_gra))*ez)*ex-sin(ang_gra)*ey;
   coe[2][1]=((1-cos(ang_gra))*ez)*ey+sin(ang_gra)*ex;
   coe[2][2]=((1-cos(ang_gra))*ez)*ex+cos(ang_gra);
}

Tupla4f Matriz4x4f::operator*( Tupla4f & tupla )const
{
   Tupla4f res;
   for( unsigned fila = 0 ; fila < 4 ; fila++ )
   {  res.coo[fila] = 0.0f ;
      for( unsigned colu = 0 ; colu < 4 ; colu++ )
         res.coo[fila] += coe[fila][colu]*tupla.coo[colu] ;
   }
   return res ;
}

Matriz4x4f Matriz4x4f::operator*( const Matriz4x4f & der )
{
   Matriz4x4f res ;
   for( unsigned fil = 0 ; fil < 4 ; fil++ )
      for( unsigned col = 0 ; col < 4 ; col++ )
      { res.coe[col][fil] = 0.0f ;
         for( unsigned k = 0 ; k < 4 ; k++ )
            res.coe[col][fil] += coe[k][fil]*der.coe[col][k];
      }
   return res ;
}

Matriz4x4f Matriz4x4f::operator=(const Matriz4x4f &matriz){
   if(&matriz!=this){
      for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
            this->coe[i][j]=matriz.coe[i][j];
   }
   return *this;
}
