#include <GL/glut.h>
#include "matriz.hpp"
#include "malla.hpp"
#include <iostream>



MallaTA::MallaTA()
{
   num_ver=0;
   num_tri=0;

   ver=NULL;
   tri=NULL;

   normalesCaras=NULL;
   normalesVertices=NULL;

   tipoDibuj=GL_LINE;
   modoAjedr=false;
   caras=false;
   vertices=false;
}


MallaTA::MallaTA(int numv, int numc)
{
   num_ver=numv;
   num_tri=numc;

	ver = new Tupla3f[num_ver];
	tri = new Tupla3i[num_tri];

   normalesCaras=NULL;
   normalesVertices=NULL;

   tipoDibuj=GL_LINE;
   modoAjedr=false;
   caras=false;
   vertices=false;
}

void MallaTA::inicializarMalla(int numv, int numc)
{
   num_ver = numv;
   num_tri = numc;

   ver = new Tupla3f[num_ver];
	tri = new Tupla3i[num_tri];
   normalesCaras = new Tupla3f[num_tri];
   normalesVertices = new Tupla3f[num_ver];

   tipoDibuj=GL_LINE;
   modoAjedr=false;
   caras=false;
   vertices=false;
}

Tupla3f MallaTA::P2_CalcularCentro(int i){
   Tupla3f coef((ver[tri[i].idx[0]].coo[0]+ver[tri[i].idx[1]].coo[0]+ver[tri[i].idx[2]].coo[0])/3.0, (ver[tri[i].idx[0]].coo[1]+ver[tri[i].idx[1]].coo[1]+ver[tri[i].idx[2]].coo[1])/3.0, (ver[tri[i].idx[0]].coo[2]+ver[tri[i].idx[1]].coo[2]+ver[tri[i].idx[2]].coo[2])/3.0);
   return coef;
}

void MallaTA::extraerVerTri(std::vector<float> verticesCirculo_ply, std::vector<int> carasCirculo_ply)
{
	for (unsigned i=0; i<num_tri; i++){
		tri[i].idx[0]=carasCirculo_ply[i*3];
		tri[i].idx[1]=carasCirculo_ply[i*3+1];
		tri[i].idx[2]=carasCirculo_ply[i*3+2];
	}

	for (unsigned i=0; i<num_ver; i++){
		ver[i].coo[0]=verticesCirculo_ply[i*3];
		ver[i].coo[1]=verticesCirculo_ply[i*3+1];
		ver[i].coo[2]=verticesCirculo_ply[i*3+2];
	}
}

void MallaTA::copiarTriangulos(MallaTA malla)
{
	for (unsigned i=0; i<num_tri; i++){
		tri[i].idx[0]=malla.tri[i].idx[0];
		tri[i].idx[1]=malla.tri[i].idx[1];
		tri[i].idx[2]=malla.tri[i].idx[2];
	}
}

void MallaTA::aplicarMatriz(Matriz4x4f matriz, MallaTA malla)
{
   for (unsigned i=0; i<num_ver; i++)
   {
      Tupla4f coor(malla.ver[i].coo[0],malla.ver[i].coo[1],malla.ver[i].coo[2],1);
      Tupla4f coord=matriz * coor;

      ver[i].coo[0]=coord.coo[0];
      ver[i].coo[1]=coord.coo[1];
      ver[i].coo[2]=coord.coo[2];
   }
}

void MallaTA::dibujarObjetos() 
{
   glPolygonMode(GL_FRONT_AND_BACK,tipoDibuj);

   glEnableClientState( GL_VERTEX_ARRAY );

   glVertexPointer( 3, GL_FLOAT, 0, ver );

   if (modoAjedr == true){
      for (unsigned i=0; i<num_tri; i++){
         if (i%2==0) glColor3f(1.0,0.0,0.0);
	else glColor3f(1.0,1.0,1.0);

	 glDrawElements( GL_TRIANGLES, 3, GL_UNSIGNED_INT, tri[i].idx ) ;
      }
   }
   else glDrawElements( GL_TRIANGLES, 3*num_tri, GL_UNSIGNED_INT, tri ) ;
   glColor3f(0.0,0.0,1.0);
   if (caras){
      for (unsigned i=0; i<num_tri; i++){
         Tupla3f coef(P2_CalcularCentro(i));
         glBegin(GL_LINES);
         glVertex3f(coef.coo[0], coef.coo[1], coef.coo[2]);
         glVertex3f(normalesCaras[i].coo[0]/7+coef.coo[0], normalesCaras[i].coo[1]/7+coef.coo[1], normalesCaras[i].coo[2]/7+coef.coo[2]);
         glEnd();
      }
   }
   glColor3f(0.0,1.0,0.0);
   if (vertices){
      for (unsigned i=0;i<num_ver;i++){
            glBegin(GL_LINES);
            glVertex3f(ver[i].coo[0], ver[i].coo[1], ver[i].coo[2]);
            glVertex3f(normalesVertices[i].coo[0]/8+ver[i].coo[0], normalesCaras[i].coo[1]/8+ver[i].coo[1], normalesCaras[i].coo[2]/8+ver[i].coo[2]);
            glEnd();
      }
   }
}

void MallaTA::cambiarCaracsDibujo(GLint tipo, bool ajedrez){
	tipoDibuj=tipo;
	modoAjedr=ajedrez;
}


// Practica 2
/*
void MallaTA::normalesCaras(bool c){
   caras=c;
}

void MallaTA::normalesVertices(bool v){
   vertices=v;
}
*/
void MallaTA::P2_CrearTriangulo(int &num_tri, int ver1, int ver2, int ver3){
   tri[num_tri].idx[0]=ver1;
   tri[num_tri].idx[1]=ver2;
   tri[num_tri].idx[2]=ver3;

   num_tri++;
}

void MallaTA::calcularNormales(){

   for (unsigned i=0; i<num_ver;i++) normalesVertices[i]=Tupla3f(0,0,0);

   for(unsigned i=0; i<num_tri;i++){
      normalesCaras[i]=(ver[tri[i].idx[1]]-ver[tri[i].idx[0]])*(ver[tri[i].idx[2]]-ver[tri[i].idx[0]]);
      normalesCaras[i]=normalesCaras[i]/sqrt(normalesCaras[i]|normalesCaras[i]);

      normalesVertices[tri[i].idx[0]]=normalesCaras[i]+normalesVertices[tri[i].idx[0]];
      normalesVertices[tri[i].idx[1]]=normalesCaras[i]+normalesVertices[tri[i].idx[1]];
      normalesVertices[tri[i].idx[2]]=normalesCaras[i]+normalesVertices[tri[i].idx[2]];
   }

   for (unsigned i=0; i<num_ver;i++) normalesVertices[i]=normalesVertices[i]/sqrt(normalesVertices[i]|normalesVertices[i]);
}


int MallaTA::P2_PreparaPerfil(Tupla3f &ver1, Tupla3f &ver2, std::vector<float> pPerfil_ply){  //Devuelve el número de puntos del perfil sin contar los extremos
   unsigned num_puntos=pPerfil_ply.size()/3;
   int sum=0;
      if (pPerfil_ply[0]==0) {
         num_puntos--;
         sum++;
         ver1=Tupla3f(pPerfil_ply[0], pPerfil_ply[1], pPerfil_ply[2]);
      }
      else ver1=Tupla3f(0, pPerfil_ply[1], pPerfil_ply[2]);
      if (pPerfil_ply[pPerfil_ply.size()-3]==0) {
         num_puntos--;
         ver2=Tupla3f(pPerfil_ply[pPerfil_ply.size()-3], pPerfil_ply[pPerfil_ply.size()-2], pPerfil_ply[pPerfil_ply.size()-1]); 
      }
      else ver2=Tupla3f(0, pPerfil_ply[pPerfil_ply.size()-2], pPerfil_ply[pPerfil_ply.size()-1]);

      for (unsigned i=0; i<num_puntos; i++){
         ver[i].coo[0]=pPerfil_ply[i*3+sum];
         ver[i].coo[1]=pPerfil_ply[i*3+1+sum];
         ver[i].coo[2]=pPerfil_ply[i*3+2+sum];
      }
      return num_puntos;
}


void MallaTA::P2_AnadirExtremos(int ultimo, int num_p, Tupla3f ver1, Tupla3f ver2){  //Incluye en el vector de vertices los extremos de la silueta
   ver[ultimo].coo[0]=ver1.coo[0];
   ver[ultimo].coo[1]=ver1.coo[1];
   ver[ultimo].coo[2]=ver1.coo[2];

   ver[ultimo+1].coo[0]=ver2.coo[0];
   ver[ultimo+1].coo[1]=ver2.coo[1];
   ver[ultimo+1].coo[2]=ver2.coo[2];
}

void MallaTA::P2_RotarPunto(int anterior1, int ultimo, const Matriz4x4f matrizRotacion){
   Tupla4f coor(ver[anterior1].coo[0],ver[anterior1].coo[1],ver[anterior1].coo[2],1);
   Tupla4f Ult=matrizRotacion*coor;

   ver[ultimo].coo[0]=Ult.coo[0];
   ver[ultimo].coo[1]=Ult.coo[1];
   ver[ultimo].coo[2]=Ult.coo[2];
}


void MallaTA::crearRotacion(std::vector<float> pPerfil_ply, int M)
{
   N=M;
   // Se crea la matriz de rotación con respecto al eje Y
   Matriz4x4f matrizRotacion;
   matrizRotacion.MAT_Rotacion( 2*M_PI/ N, 0, 1, 0 );

   Tupla3f ver1;
   Tupla3f ver2;

   unsigned num_p=P2_PreparaPerfil(ver1,ver2, pPerfil_ply);

   int anterior1=0, anterior2;
   int ultimo=num_p;
   int num_triangulo=0;

   for (unsigned j=0; j<N-1; j++){
      for (unsigned i=0; i<num_p-2; i++){

         for (int k=0; k<2;k++){
            if (i==0 || k==0) P2_RotarPunto(anterior1, ultimo, matrizRotacion);
            if (i==0 && k==0){
               P2_CrearTriangulo(num_triangulo, anterior1, anterior1+1, ultimo);
               anterior2=ultimo;
               ultimo++;
               anterior1++;
            }
         }

         P2_CrearTriangulo(num_triangulo, anterior1, ultimo, anterior2);
         P2_CrearTriangulo(num_triangulo, anterior1, anterior1+1, ultimo);
         anterior1++;
         anterior2=ultimo;
         ultimo++;

         if (i==num_p-3){
            P2_RotarPunto(anterior1, ultimo, matrizRotacion);
            P2_CrearTriangulo(num_triangulo, anterior1, ultimo, anterior2);
            anterior1++;
            ultimo++;
         }
      }
   }


   anterior2=0;
   for (unsigned i=0; i<num_p-1; i++){
      P2_CrearTriangulo(num_triangulo, anterior2, anterior1, anterior1+1);
      P2_CrearTriangulo(num_triangulo, anterior2, anterior1+1, anterior2+1);
      anterior1++;
      anterior2++;
   }

   P2_AnadirExtremos(ultimo,num_p, ver1, ver2);

   for (unsigned i=0;i<N;i++){
      if (i!=N-1) tri[num_triangulo].idx[2]= num_p*(i+1);
      else tri[num_triangulo].idx[2]=0;

      P2_CrearTriangulo(num_triangulo, ultimo, num_p*i, tri[num_triangulo].idx[2]);
   }
   for (unsigned i=0;i<N;i++){
      if (i!=N-1) tri[num_triangulo].idx[2]= num_p*(i+1)+num_p-1;
      else tri[num_triangulo].idx[2]=num_p-1;

      P2_CrearTriangulo(num_triangulo, ultimo+1, tri[num_triangulo].idx[2], num_p*i+num_p-1);
   }
}

