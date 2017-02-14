//
//
//++++++++++++++++++++++++++++nonlinear interpolation for all models++++++++++++++++++++++++++++++++++++++
//
//

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <vector>
#include "objMesh.h"//
#include <gl/freeglut.h>
#include "mat3d.h"
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"//

using namespace Eigen;
using namespace std;

//===========Number of vertices/ faces=======

unsigned int NumV = 0;
unsigned int NumF = 0;
unsigned int NumG = 0; 
int NumUnE=0;
double TimeStep=0.25;


//====================Example-space==================

ObjMesh *MESH;
ObjMesh *example_1;
ObjMesh *example_2;


//===================================================Half-edge===================================================
struct HE_Edge;
struct HE_Face;
struct HE_Vert;

struct HE_Edge
{
	HE_Vert *startvert;
	HE_Face *leftface;
	HE_Edge *pair;
	HE_Edge *pre;
	HE_Edge *next;
	double dihedral_angle;
	double length;
	double area_sum_neighbor;
	int eindex;
};	

struct HE_Face
{
	HE_Edge *firstedge;
	int faceindex;
	Vec3d facenorm;
	double facearea;
};	

struct HE_Vert
{
	Vec3d VertPosition;
	HE_Edge *outedge;
	int vertindex;
	Vec3d VertNorm;
	bool isCon;
};

struct Edge_list
{
	int HEedgeID;
	bool isEdge;
};

vector<HE_Edge> HEedge;
vector<HE_Face> HEface;
vector<HE_Vert> HEvert;
HE_Edge ***isPair;//original

vector<HE_Face> HEface_curr;
vector<HE_Vert> HEvert_curr;
vector<int> HEvert_Con;
vector<double> Dihedral_curr;//position after move

vector<HE_Face> HEface_e1;
vector<HE_Vert> HEvert_e1;
vector<double> Dihedral_e1;//position of example_1

vector<HE_Face> HEface_e2;
vector<HE_Vert> HEvert_e2;
vector<double> Dihedral_e2;//position of example_1

vector<vector<Edge_list> > Edgelist;
vector<int> hedgeid;

int MultiMode=0;

class Vertex
{
public:
	Vec3d vertpos;
	bool selected;
	bool iscon;
	
	Vertex()
		: vertpos(0,0,0)
		, selected(false),iscon(false)
	{}

	void Draw()
	{	
		if(selected)
		{
			if (MultiMode>0)
			{
				glColor3f(1.0,0.0,1.0);
				glPushMatrix();
					glTranslated(vertpos.operator[](0),vertpos.operator[](1),vertpos.operator[](2));
					glutSolidSphere(0.05,10,10);
				glPopMatrix();//
			}
			else
			{
				glColor3f(0.0,1.0,0.0);
				glPushMatrix();
					glTranslated(vertpos.operator[](0),vertpos.operator[](1),vertpos.operator[](2));
					glutSolidSphere(0.05,10,10);
				glPopMatrix();//
			}
		
		}
		else if (iscon)
		{
			glColor3f(1.0,1.0,1.0);
		}
		else
		{
			glColor3f(1.0,0.0,0.0);
		}
		glPointSize(10.0f);
		glBegin(GL_POINTS);	
			glVertex3d(vertpos.operator[](0),vertpos.operator[](1),vertpos.operator[](2));
		glEnd();
	}
};
vector<Vertex> vert_picking;

void HEvertCurr()
{
	HEvert_curr.resize(NumV);
	for(unsigned int iGroup=0;iGroup<NumG;iGroup++)
	{
		for(unsigned int iFace=0;iFace<NumF;iFace++)
		{
			Vec3d vec1_curr=Vec3d(0,0,0);
			Vec3d vec2_curr=Vec3d(0,0,0);
			Vec3d vec3_curr=Vec3d(0,0,0);
			Vec3d norm_temp=Vec3d(0,0,0);
	
			unsigned int vid0=MESH->getVertexIndex(iGroup,iFace,0);
			unsigned int vid1=MESH->getVertexIndex(iGroup,iFace,1);
			unsigned int vid2=MESH->getVertexIndex(iGroup,iFace,2);//assume all faces are triangle
	
			HEvert_curr[vid0].vertindex=HEvert[vid0].vertindex;
			HEvert_curr[vid0].VertPosition=HEvert[vid0].VertPosition+TimeStep*(HEvert_e1[vid0].VertPosition-HEvert[vid0].VertPosition);
			HEvert_curr[vid0].isCon=false;
			//HEvert_curr[vid0].outedge=HEvert[vid0].outedge;//
			HEvert_curr[vid1].vertindex=HEvert[vid1].vertindex;
			HEvert_curr[vid1].VertPosition=HEvert[vid1].VertPosition+TimeStep*(HEvert_e1[vid1].VertPosition-HEvert[vid1].VertPosition);
			HEvert_curr[vid1].isCon=false;
			//HEvert_curr[vid1].outedge=HEvert[vid1].outedge;//
			HEvert_curr[vid2].vertindex=HEvert[vid2].vertindex;
			HEvert_curr[vid2].VertPosition=HEvert[vid2].VertPosition+TimeStep*(HEvert_e1[vid2].VertPosition-HEvert[vid2].VertPosition);
			HEvert_curr[vid2].isCon=false;
		}
	}

	for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
	{
		HEvert_curr[HEvert_Con[ivert]].isCon=true;
		HEvert[HEvert_Con[ivert]].isCon=true;
		HEvert_curr[HEvert_Con[ivert]].VertPosition=HEvert[HEvert_Con[ivert]].VertPosition;
	}
}

void HEfaceCurr()
{
	HEface_curr.resize(NumF);
	for(unsigned int iGroup=0;iGroup<NumG;iGroup++)
	{
		for(unsigned int iFace=0;iFace<NumF;iFace++)
		{
			Vec3d vec1_curr=Vec3d(0,0,0);
			Vec3d vec2_curr=Vec3d(0,0,0);
			Vec3d vec3_curr=Vec3d(0,0,0);
			Vec3d norm_temp=Vec3d(0,0,0);

			unsigned int vid0=MESH->getVertexIndex(iGroup,iFace,0);
			unsigned int vid1=MESH->getVertexIndex(iGroup,iFace,1);
			unsigned int vid2=MESH->getVertexIndex(iGroup,iFace,2);//assume all faces are triangle

			vec1_curr=HEvert_curr[vid1].VertPosition-HEvert_curr[vid0].VertPosition;
			vec2_curr=HEvert_curr[vid2].VertPosition-HEvert_curr[vid0].VertPosition;
			vec3_curr=HEvert_curr[vid2].VertPosition-HEvert_curr[vid1].VertPosition;
			norm_temp=cross(vec1_curr,vec3_curr); //
			HEface_curr[iFace].faceindex=iFace;
			HEface_curr[iFace].facenorm=norm(norm_temp);
			HEface_curr[iFace].facearea=0.5*sqrt(dot(norm_temp,norm_temp));//
		}
	}
}

void DihedralCurr()
{
	Dihedral_curr.resize(3*NumF);
	for(unsigned int iedge=0;iedge<3*NumF;iedge++)
	{
		double DOT=dot((HEface_curr[HEedge[iedge].leftface->faceindex].facenorm),(HEface_curr[HEedge[iedge].pair->leftface->faceindex].facenorm));
		if (abs(DOT)<1)
		{
			Dihedral_curr[iedge]=M_PI-acos(DOT);
		}
		else
		{
			Dihedral_curr[iedge]=M_PI;
		}
	}
}

void HEvnormCurr()
{
	for(unsigned int ivert=0;ivert<NumV;ivert++)//computation of dihedral angle and neighboring area with respect to each edge
	{
		HE_Edge *out=HEvert[ivert].outedge;
		HE_Edge *curr=HEvert[ivert].outedge;
		Vec3d temp_facenorm=Vec3d(0,0,0);
		double area=0;//
		do
		{
			temp_facenorm+=HEface_curr[curr->leftface->faceindex].facenorm;
			area+=HEface_curr[curr->leftface->faceindex].facearea;
			curr=curr->pair->next;
		}while(curr!=out);
		HEvert_curr[ivert].VertNorm=norm(temp_facenorm/area);//
	}
}

void HalfEdge_MESH()
{
	int HEedgeNum=3*NumF;
	int iEdge=0;

	HEedge.resize(HEedgeNum);
	HEface.resize(NumF);
	HEvert.resize(NumV);
	
	isPair = new HE_Edge**[NumV];
	for (int i= 0; i<NumV;i++)
	{
		isPair[i] = new HE_Edge*[20];  //assume the edges connect to one point is not more than 20
	}
	for (int i= 0; i<NumV;i++)
	{
		for (int j= 0; j<20;j++)
		{
			isPair[i][j]= NULL;
		}
	}

	for(unsigned int iGroup=0;iGroup<NumG;iGroup++)
	{
		for(unsigned int iFace=0;iFace<NumF;iFace++)
		{
			Vec3d vec1=Vec3d(0,0,0);
			Vec3d vec2=Vec3d(0,0,0);
			Vec3d vec3=Vec3d(0,0,0);
			Vec3d vec1_curr=Vec3d(0,0,0);
			Vec3d vec2_curr=Vec3d(0,0,0);
			Vec3d vec3_curr=Vec3d(0,0,0);
			double length01=0;
			double length02=0;
			double length12=0;
			Vec3d norm_temp=Vec3d(0,0,0);
		
			unsigned int vid0=MESH->getVertexIndex(iGroup,iFace,0);
			unsigned int vid1=MESH->getVertexIndex(iGroup,iFace,1);
			unsigned int vid2=MESH->getVertexIndex(iGroup,iFace,2);//assume all faces are triangle

			HEvert[vid0].vertindex=vid0;
			HEvert[vid0].VertPosition=MESH->getPosition(vid0);
			HEvert[vid0].isCon=false;
			HEvert[vid1].vertindex=vid1;
			HEvert[vid1].VertPosition=MESH->getPosition(vid1);
			HEvert[vid1].isCon=false;
			HEvert[vid2].vertindex=vid2;
			HEvert[vid2].VertPosition=MESH->getPosition(vid2);
			HEvert[vid2].isCon=false;//vert

			vec1=HEvert[vid1].VertPosition-HEvert[vid0].VertPosition;
			//vec2=HEvert[vid2].VertPosition-HEvert[vid0].VertPosition;
			vec3=HEvert[vid2].VertPosition-HEvert[vid1].VertPosition;
			norm_temp=cross(vec1,vec3); 
			HEface[iFace].faceindex=iFace;
			HEface[iFace].firstedge=&HEedge[iEdge];
			HEface[iFace].facenorm=norm(norm_temp);
			HEface[iFace].facearea=0.5*sqrt(dot(norm_temp,norm_temp));//face

			HEedge[iEdge].leftface=&HEface[iFace];
			HEedge[iEdge].startvert=&HEvert[vid0];
			HEedge[iEdge].eindex=iEdge;
			//HEedge[iEdge].pair=NULL;
			HEvert[vid0].outedge=&HEedge[iEdge];
			for (int p = 0; p < 20; p++)
			{
				if (isPair[vid1][p]==NULL)
				{
					isPair[vid1][p]=&HEedge[iEdge];
					break;
				}
			}
			for (int q = 0; q < 20; q++)
			{
				if (isPair[vid0][q] != NULL)
				{
					if (isPair[vid0][q]->startvert ==&HEvert[vid1] )
					{
						isPair[vid0][q]->pair = &HEedge[iEdge];
						HEedge[iEdge].pair = isPair[vid0][q];
						break;
					}
				}
				else
				{
					break;
				}
			}///

			HEedge[iEdge+1].leftface=&HEface[iFace];
			HEedge[iEdge+1].startvert=&HEvert[vid1];
			HEedge[iEdge+1].eindex=iEdge+1;
			//HEedge[iEdge+1].pair=NULL;
			HEvert[vid1].outedge=&HEedge[iEdge+1];
			for (int p = 0; p < 20; p++)
			{
				if (isPair[vid2][p]==NULL)
				{
					isPair[vid2][p]=&HEedge[iEdge+1];
					break;
				}
			}
			for (int q = 0; q < 20; q++)
			{
				if (isPair[vid1][q] != NULL)
				{
					if (isPair[vid1][q]->startvert ==&HEvert[vid2] )
					{
						isPair[vid1][q]->pair = &HEedge[iEdge+1];
						HEedge[iEdge+1].pair = isPair[vid1][q];
						break;
					}
				}
				else
				{
					break;
				}
			}///

			HEedge[iEdge+2].leftface=&HEface[iFace];
			HEedge[iEdge+2].startvert=&HEvert[vid2];
			HEedge[iEdge+2].eindex=iEdge+2;
			//HEedge[iEdge+2].pair=NULL;
			HEvert[vid2].outedge=&HEedge[iEdge+2];
			for (int p = 0; p < 20; p++)
			{
				if (isPair[vid0][p]==NULL)
				{
					isPair[vid0][p]=&HEedge[iEdge+2];
					break;
				}
			}
			for (int q = 0; q < 20; q++)
			{
				if (isPair[vid2][q] != NULL)
				{
					if (isPair[vid2][q]->startvert ==&HEvert[vid0] )
					{
						isPair[vid2][q]->pair = &HEedge[iEdge+2];
						HEedge[iEdge+2].pair = isPair[vid2][q];
						break;
					}
				}
				else
				{
					break;
				}
			}///

			HEedge[iEdge].pre=&HEedge[iEdge+2];
			HEedge[iEdge].next=&HEedge[iEdge+1];
			HEedge[iEdge+1].pre=&HEedge[iEdge];
			HEedge[iEdge+1].next=&HEedge[iEdge+2];
			HEedge[iEdge+2].pre=&HEedge[iEdge+1];
			HEedge[iEdge+2].next=&HEedge[iEdge];//
			
			Vec3d vec_length;
			vec_length=(HEedge[iEdge].next->startvert->VertPosition)-(HEedge[iEdge].startvert->VertPosition);
			HEedge[iEdge].length=len(vec_length);
			vec_length=(HEedge[iEdge+1].next->startvert->VertPosition)-(HEedge[iEdge+1].startvert->VertPosition);
			HEedge[iEdge+1].length=len(vec_length);
			vec_length=(HEedge[iEdge+2].next->startvert->VertPosition)-(HEedge[iEdge+2].startvert->VertPosition);
			HEedge[iEdge+2].length=len(vec_length);//edge length

			//store the edgelist
			if ((HEedge[iEdge].startvert->vertindex)<(HEedge[iEdge].next->startvert->vertindex))
			{
				Edgelist[HEedge[iEdge].startvert->vertindex][HEedge[iEdge].next->startvert->vertindex].HEedgeID=iEdge;
				//Edgelist[HEedge[iEdge].startvert->vertindex][HEedge[iEdge].next->startvert->vertindex].length=HEedge[iEdge].length;
				Edgelist[HEedge[iEdge].startvert->vertindex][HEedge[iEdge].next->startvert->vertindex].isEdge=true;
			}
			if ((HEedge[iEdge+1].startvert->vertindex)<(HEedge[iEdge+1].next->startvert->vertindex))
			{
				Edgelist[HEedge[iEdge+1].startvert->vertindex][HEedge[iEdge+1].next->startvert->vertindex].HEedgeID=iEdge+1;
				//Edgelist[HEedge[iEdge+1].startvert->vertindex][HEedge[iEdge+1].next->startvert->vertindex].length=HEedge[iEdge+1].length;
				Edgelist[HEedge[iEdge+1].startvert->vertindex][HEedge[iEdge+1].next->startvert->vertindex].isEdge=true;
			}
			if ((HEedge[iEdge+2].startvert->vertindex)<(HEedge[iEdge+2].next->startvert->vertindex))
			{
				Edgelist[HEedge[iEdge+2].startvert->vertindex][HEedge[iEdge+2].next->startvert->vertindex].HEedgeID=iEdge+2;
				//Edgelist[HEedge[iEdge+2].startvert->vertindex][HEedge[iEdge+2].next->startvert->vertindex].length=HEedge[iEdge+2].length;
				Edgelist[HEedge[iEdge+2].startvert->vertindex][HEedge[iEdge+2].next->startvert->vertindex].isEdge=true;
			}

			iEdge=iEdge+3;
		}//////end for-loop of faces
		
		for(unsigned int iedge=0;iedge<3*NumF;iedge++)//dihedral angle and neighboring area
		{
			double DOT=dot((HEedge[iedge].leftface->facenorm),(HEedge[iedge].pair->leftface->facenorm));
			if (abs(DOT)<1)
			{
				HEedge[iedge].dihedral_angle=M_PI-acos(DOT);
			}
			else
			{
				HEedge[iedge].dihedral_angle=M_PI;
			}
			HEedge[iedge].area_sum_neighbor=HEedge[iedge].leftface->facearea+HEedge[iedge].pair->leftface->facearea;
		}//////loop for dihedral angle and neighboring area

		for(unsigned int ivert=0;ivert<NumV;ivert++)//vert normal
		{
			HE_Edge *out=HEvert[ivert].outedge;
			HE_Edge *curr=HEvert[ivert].outedge;
			Vec3d temp_facenorm=Vec3d(0,0,0);
			double area=0;//
			do
			{
				temp_facenorm+=HEface[curr->leftface->faceindex].facenorm;
				area+=HEface[curr->leftface->faceindex].facearea;
				curr=curr->pair->next;
			}while(curr!=out);
			HEvert[ivert].VertNorm=norm(temp_facenorm/area);//
		}//loop for vert normal
	}//end for-loop of groups
}//

void HalfEdge_example_1()  //vert and face
{
	Vec3d vec1=Vec3d(0,0,0);
	Vec3d vec2=Vec3d(0,0,0);

	HEface_e1.resize(NumF);
	HEvert_e1.resize(NumV);
	Dihedral_e1.resize(3*NumF);

	for(unsigned int iGroup=0;iGroup<NumG;iGroup++)
	{
		for(unsigned int iFace=0;iFace<NumF;iFace++)
		{
			Vec3d vec1=Vec3d(0,0,0);
			Vec3d vec2=Vec3d(0,0,0);
			Vec3d vec3=Vec3d(0,0,0);
			double length01=0;
			double length02=0;
			double length12=0;
			Vec3d norm_temp=Vec3d(0,0,0);
		
			unsigned int vid0=example_1->getVertexIndex(iGroup,iFace,0);
			unsigned int vid1=example_1->getVertexIndex(iGroup,iFace,1);
			unsigned int vid2=example_1->getVertexIndex(iGroup,iFace,2);//assume all faces are triangle

			HEvert_e1[vid0].vertindex=vid0;
			HEvert_e1[vid0].VertPosition=example_1->getPosition(vid0);
			HEvert_e1[vid1].vertindex=vid1;
			HEvert_e1[vid1].VertPosition=example_1->getPosition(vid1);
			HEvert_e1[vid2].vertindex=vid2;
			HEvert_e1[vid2].VertPosition=example_1->getPosition(vid2);//vert

			vec1=HEvert_e1[vid1].VertPosition-HEvert_e1[vid0].VertPosition;
			vec2=HEvert_e1[vid2].VertPosition-HEvert_e1[vid0].VertPosition;
			vec3=HEvert_e1[vid2].VertPosition-HEvert_e1[vid1].VertPosition;
			norm_temp=cross(vec1,vec3); 
			HEface_e1[iFace].faceindex=iFace;
			//HEface_e1[iFace].firstedge=&HEedge_e1[iEdge];
			HEface_e1[iFace].facenorm=norm(norm_temp);
			HEface_e1[iFace].facearea=0.5*sqrt(dot(norm_temp,norm_temp));//face

		}//end for-loop of faces
		
		for(unsigned int iedge=0;iedge<3*NumF;iedge++)//computation of dihedral angle and neighboring area with respect to each edge
		{
			unsigned int faceid1=HEedge[iedge].leftface->faceindex;
			unsigned int faceid2=HEedge[iedge].pair->leftface->faceindex;
			double DOT=dot((HEface_e1[faceid1].facenorm),(HEface_e1[faceid2].facenorm));
			if (abs(DOT)<1)
			{
				Dihedral_e1[iedge]=M_PI-acos(DOT);
			}
			else
			{
				Dihedral_e1[iedge]=M_PI;
			}
		}

		for(unsigned int ivert=0;ivert<NumV;ivert++)//vert normal
		{
			HE_Edge *out=HEvert[ivert].outedge;
			HE_Edge *curr=HEvert[ivert].outedge;
			Vec3d temp_facenorm=Vec3d(0,0,0);
			double area=0;//
			do
			{
				temp_facenorm+=HEface_e1[curr->leftface->faceindex].facenorm;
				area+=HEface_e1[curr->leftface->faceindex].facearea;
				curr=curr->pair->next;
			}while(curr!=out);
			HEvert_e1[ivert].VertNorm=norm(temp_facenorm/area);//
		}//loop for vert normal

	}//end for-loop of groups
}//

void HalfEdge_example_2()  //vert and face
{
	Vec3d vec1=Vec3d(0,0,0);
	Vec3d vec2=Vec3d(0,0,0);

	HEface_e2.resize(NumF);
	HEvert_e2.resize(NumV);
	Dihedral_e2.resize(3*NumF);

	for(unsigned int iGroup=0;iGroup<NumG;iGroup++)
	{
		for(unsigned int iFace=0;iFace<NumF;iFace++)
		{
			Vec3d vec1=Vec3d(0,0,0);
			Vec3d vec2=Vec3d(0,0,0);
			Vec3d vec3=Vec3d(0,0,0);
			double length01=0;
			double length02=0;
			double length12=0;
			Vec3d norm_temp=Vec3d(0,0,0);
		
			unsigned int vid0=example_2->getVertexIndex(iGroup,iFace,0);
			unsigned int vid1=example_2->getVertexIndex(iGroup,iFace,1);
			unsigned int vid2=example_2->getVertexIndex(iGroup,iFace,2);//assume all faces are triangle

			HEvert_e2[vid0].vertindex=vid0;
			HEvert_e2[vid0].VertPosition=example_2->getPosition(vid0);
			HEvert_e2[vid1].vertindex=vid1;
			HEvert_e2[vid1].VertPosition=example_2->getPosition(vid1);
			HEvert_e2[vid2].vertindex=vid2;
			HEvert_e2[vid2].VertPosition=example_2->getPosition(vid2);//vert

			vec1=HEvert_e2[vid1].VertPosition-HEvert_e2[vid0].VertPosition;
			vec2=HEvert_e2[vid2].VertPosition-HEvert_e2[vid0].VertPosition;
			vec3=HEvert_e2[vid2].VertPosition-HEvert_e2[vid1].VertPosition;
			norm_temp=cross(vec1,vec3); 
			HEface_e2[iFace].faceindex=iFace;
			//HEface_e1[iFace].firstedge=&HEedge_e1[iEdge];
			HEface_e2[iFace].facenorm=norm(norm_temp);
			HEface_e2[iFace].facearea=0.5*sqrt(dot(norm_temp,norm_temp));//face

		}//end for-loop of faces
		
		for(unsigned int iedge=0;iedge<3*NumF;iedge++)//computation of dihedral angle and neighboring area with respect to each edge
		{
			unsigned int faceid1=HEedge[iedge].leftface->faceindex;
			unsigned int faceid2=HEedge[iedge].pair->leftface->faceindex;
			Dihedral_e2[iedge]=M_PI-acos(dot((HEface_e2[faceid1].facenorm),(HEface_e2[faceid2].facenorm)));

		}
	}//end for-loop of groups
}//


//=======================================================Weight=======================================================

vector<double> weight_s,weight_b;
double weight_v;
double para_v=1000, error=0.01;// para_lamda, para_u,
vector<double> para_lamda,para_u;
double para_lamda_temp=10000, para_u_temp=1;
double volume,volume_curr,volume_e1,volume_e2;

double Volume(ObjMesh* mesh,vector<HE_Vert> &vert)
{
	double vol = 0;
	for(unsigned int iGroup=0;iGroup<NumG;iGroup++)
	{
		for(unsigned int iFace=0;iFace<NumF;iFace++)
		{
			double volume_temp=0;
			
			int vert0=mesh->getVertexIndex(iGroup,iFace,0);
			int vert1=mesh->getVertexIndex(iGroup,iFace,1);
			int vert2=mesh->getVertexIndex(iGroup,iFace,2);

			Vec3d vertposition0=vert[vert0].VertPosition;
			Vec3d vertposition1=vert[vert1].VertPosition;
			Vec3d vertposition2=vert[vert2].VertPosition;
			
			volume_temp=(abs(dot((cross(vertposition0,vertposition1)),vertposition2)))/6;
			vol+=volume_temp;
		}
	}
	return vol;
}//loop for volume 

void Weight()
{
	int i=0;
	volume = Volume(MESH,HEvert);
	weight_v=sqrt(para_v)/volume;
	//cout<<"weight_v"<<"\t"<<weight_v<<endl;//

	for (int iedge= 0; iedge < (1.5*NumF); iedge++)
	{
		if (hedgeid[iedge]>=0)
		{
			//cout<<"hedgeid"<<iedge<<"\t"<<hedgeid[iedge]<<endl;
			weight_s.emplace_back(sqrt(para_lamda_temp)/HEedge[hedgeid[iedge]].length);
			weight_b.emplace_back((sqrt(para_u_temp)*HEedge[hedgeid[iedge]].length)/(sqrt(HEedge[hedgeid[iedge]].area_sum_neighbor)));//
		}
	}
}


//===============================================Mesh Deformation======================================================

VectorXd ResiFunc;
SparseMatrix <double> Jaco_j;
VectorXd Jaco_u;
vector< Triplet<double> > triplets_j; 
double interp1=0;

void  Resi_Func()
{
	ResiFunc.resize(2*NumUnE+1);

	int iResi=0;

	for (int i = 0; i < (3*NumF);)
	{
		if (i<(1.5*NumF))
		{
			if (hedgeid[i]>=0)
			{
				//cout<<"hedgeid in ResiFunc WS"<<i<<"\t"<<hedgeid[i]<<endl;
				double length_curr=len((HEvert_curr[HEedge[hedgeid[i]].next->startvert->vertindex].VertPosition)-(HEvert_curr[HEedge[hedgeid[i]].startvert->vertindex].VertPosition));
				//double length_curr=len((vert_picking[HEedge[hedgeid[i]].next->startvert->vertindex].vertpos)-(vert_picking[HEedge[hedgeid[i]].startvert->vertindex].vertpos));
				double length=len((HEedge[hedgeid[i]].next->startvert->VertPosition)-(HEedge[hedgeid[i]].startvert->VertPosition));
				double length_e1=len((HEvert_e1[HEedge[hedgeid[i]].next->startvert->vertindex].VertPosition)-(HEvert_e1[HEedge[hedgeid[i]].startvert->vertindex].VertPosition));
				ResiFunc[iResi]=weight_s[iResi]*(length_curr-length-(interp1*(length_e1-length)));
				//cout<<"hedgeid in ResiFunc WS"<<i<<"\t"<<hedgeid[i]<<endl;
				//cout<<"ResiFunc WS"<<iResi<<"\t"<<ResiFunc[iResi]<<endl;
				
				iResi++;
			}
			i++;
		}
		else//
		{
			//cout<<"hedgeid in ResiFunc WB"<<i<<"\t"<<hedgeid[i-1.5*NumF]<<endl;
			if (hedgeid[i-1.5*NumF]>=0)
			{
				ResiFunc[iResi]=weight_b[iResi-NumUnE]*(Dihedral_curr[hedgeid[i-1.5*NumF]]-HEedge[hedgeid[i-1.5*NumF]].dihedral_angle-(interp1*(Dihedral_e1[hedgeid[i-1.5*NumF]]-HEedge[hedgeid[i-1.5*NumF]].dihedral_angle)));		
				//cout<<"hedgeid in ResiFunc WB"<<i<<"\t"<<hedgeid[i-1.5*NumF]<<endl;
				//cout<<"ResiFunc WB"<<iResi<<"\t"<<ResiFunc[iResi]<<endl;
				iResi++;
			}
			i++;
		}
	}

	volume_curr=Volume(MESH,HEvert_curr);
	//cout<<"volume_curr"<<volume_curr<<endl;
	ResiFunc[2*weight_s.size()]=weight_v*(volume_curr-volume-interp1*(volume_e1-volume));
	//cout<<"ResiFunc"<<2*weight_s.size()<<"\t"<<ResiFunc[2*weight_s.size()]<<endl;
}

void Jaco_ju()
{
	Jaco_j.resize(2*NumUnE,(3*NumV+1));
	int i=0;
	for (int iJj = 0; iJj <(1.5*NumF); iJj++)
	{
		if (hedgeid[iJj]>=0)
		{
			Vec3d x1_a=HEvert_curr[HEedge[hedgeid[iJj]].startvert->vertindex].VertPosition;
			Vec3d x2_a=HEvert_curr[HEedge[hedgeid[iJj]].next->startvert->vertindex].VertPosition;
			Vec3d x3_a=HEvert_curr[HEedge[hedgeid[iJj]].pair->pre->startvert->vertindex].VertPosition;
			Vec3d x4_a=HEvert_curr[HEedge[hedgeid[iJj]].pre->startvert->vertindex].VertPosition;
			Vec3d e=x2_a-x1_a;////

			double EdgeLength=len((HEedge[hedgeid[iJj]].next->startvert->VertPosition)-(HEedge[hedgeid[iJj]].startvert->VertPosition));
			double EdgeLength_e1=len((HEvert_e1[HEedge[hedgeid[iJj]].next->startvert->vertindex].VertPosition)-(HEvert_e1[HEedge[hedgeid[iJj]].startvert->vertindex].VertPosition));////
			int col1_edge=3*((HEedge[hedgeid[iJj]].startvert->vertindex));
			int col2_edge=3*((HEedge[hedgeid[iJj]].next->startvert->vertindex));

			if (!(HEedge[hedgeid[iJj]].startvert->isCon))
			{
				//cout<<"Jaco_j "<<i<<"vertIndex"<<"\t"<<HEedge[hedgeid[iJj]].startvert->vertindex<<endl;
				Vec3d Vec1_edge=(-1)*weight_s[i]*(norm(e));
				triplets_j.emplace_back(i,col1_edge, Vec1_edge.operator[](0));
				triplets_j.emplace_back(i,col1_edge+1, Vec1_edge.operator[](1));
				triplets_j.emplace_back(i,col1_edge+2, Vec1_edge.operator[](2));//
			}
			if (!(HEedge[hedgeid[iJj]].next->startvert->isCon))
			{
				//cout<<"Jaco_j "<<i<<"vertIndex"<<"\t"<<HEedge[hedgeid[iJj]].next->startvert->vertindex<<endl;
				Vec3d Vec2_edge=weight_s[i]*(norm(e));
				triplets_j.emplace_back(i,col2_edge, Vec2_edge.operator[](0));
				triplets_j.emplace_back(i,col2_edge+1, Vec2_edge.operator[](1));
				triplets_j.emplace_back(i,col2_edge+2, Vec2_edge.operator[](2));//
			}
			triplets_j.emplace_back(i,(3*NumV),((-1)*weight_s[i]*(EdgeLength_e1-EdgeLength)));////edge length, weight_s

			Vec3d n1=HEface_curr[HEedge[hedgeid[iJj]].leftface->faceindex].facenorm;
			Vec3d n2=HEface_curr[HEedge[hedgeid[iJj]].pair->leftface->faceindex].facenorm;////????????????????注意检查face_curr是否在drag后更新
			int col1_a=3*(HEedge[hedgeid[iJj]].startvert->vertindex);
			int col2_a=3*(HEedge[hedgeid[iJj]].next->startvert->vertindex);
			int col3_a=3*(HEedge[hedgeid[iJj]].pair->pre->startvert->vertindex);
			int col4_a=3*(HEedge[hedgeid[iJj]].pre->startvert->vertindex);
			double xe_32,xe_31,xe_41,xe_42;
			xe_32=((x3_a.operator[](0)-x2_a.operator[](0))*(x2_a.operator[](0)-x1_a.operator[](0)))
				  +((x3_a.operator[](1)-x2_a.operator[](1))*(x2_a.operator[](1)-x1_a.operator[](1)))
				  +((x3_a.operator[](2)-x2_a.operator[](2))*(x2_a.operator[](2)-x1_a.operator[](2)));//
			xe_31=((x3_a.operator[](0)-x1_a.operator[](0))*(x2_a.operator[](0)-x1_a.operator[](0)))
				  +((x3_a.operator[](1)-x1_a.operator[](1))*(x2_a.operator[](1)-x1_a.operator[](1)))
				  +((x3_a.operator[](2)-x1_a.operator[](2))*(x2_a.operator[](2)-x1_a.operator[](2)));//
			xe_41=((x4_a.operator[](0)-x1_a.operator[](0))*(x2_a.operator[](0)-x1_a.operator[](0)))
				  +((x4_a.operator[](1)-x1_a.operator[](1))*(x2_a.operator[](1)-x1_a.operator[](1)))
				  +((x4_a.operator[](2)-x1_a.operator[](2))*(x2_a.operator[](2)-x1_a.operator[](2)));//
			xe_42=((x4_a.operator[](0)-x2_a.operator[](0))*(x2_a.operator[](0)-x1_a.operator[](0)))
				  +((x4_a.operator[](1)-x2_a.operator[](1))*(x2_a.operator[](1)-x1_a.operator[](1)))
				  +((x4_a.operator[](2)-x2_a.operator[](2))*(x2_a.operator[](2)-x1_a.operator[](2)));//????????????????实在不行检查这里是否有错
			if (!(HEedge[hedgeid[iJj]].startvert->isCon))
			{
				//cout<<"Jaco_j "<<i+NumUnE<<"vertIndex"<<"\t"<<HEedge[hedgeid[iJj]].startvert->vertindex<<endl;
				Vec3d Vec1_a=(weight_b[i]*xe_32*n1/(len(e)*len2(n1)))+(weight_b[i]*xe_42*n2/(len(e)*len2(n2)));
				triplets_j.emplace_back(i+NumUnE,col1_a, Vec1_a.operator[](0));
				triplets_j.emplace_back(i+NumUnE,col1_a+1, Vec1_a.operator[](1));
				triplets_j.emplace_back(i+NumUnE,col1_a+2, Vec1_a.operator[](2));//
			}
			if (!(HEedge[hedgeid[iJj]].next->startvert))
			{
				//cout<<"Jaco_j "<<i+NumUnE<<"vertIndex"<<"\t"<<HEedge[hedgeid[iJj]].next->startvert->vertindex<<endl;
				Vec3d Vec2_a=(-weight_b[i]*xe_31*n1/(len(e)*len2(n1)))-(weight_b[i]*xe_41*n2/(len(e)*len2(n2)));
				triplets_j.emplace_back(i+NumUnE,col2_a, Vec2_a.operator[](0));
				triplets_j.emplace_back(i+NumUnE,col2_a+1, Vec2_a.operator[](1));
				triplets_j.emplace_back(i+NumUnE,col2_a+2, Vec2_a.operator[](2));//
			}
			if (!(HEedge[hedgeid[iJj]].pair->pre->startvert->isCon))
			{
				//cout<<"Jaco_j "<<i+NumUnE<<"vertIndex"<<"\t"<<HEedge[hedgeid[iJj]].pair->pre->startvert->vertindex<<endl;
				Vec3d Vec3_a=weight_b[i]*n2*(len(e)/len2(n2));
				triplets_j.emplace_back(i+NumUnE,col3_a, Vec3_a.operator[](0));
				triplets_j.emplace_back(i+NumUnE,col3_a+1, Vec3_a.operator[](1));
				triplets_j.emplace_back(i+NumUnE,col3_a+2, Vec3_a.operator[](2));//
			}
			if (!(HEedge[hedgeid[iJj]].pre->startvert->isCon))
			{
				//cout<<"Jaco_j "<<i+NumUnE<<"vertIndex"<<"\t"<<HEedge[hedgeid[iJj]].pre->startvert->vertindex<<endl;
				Vec3d Vec4_a=weight_b[i]*n1*(len(e)/len2(n1));
				triplets_j.emplace_back(i+NumUnE,col4_a, Vec4_a.operator[](0));
				triplets_j.emplace_back(i+NumUnE,col4_a+1, Vec4_a.operator[](1));
				triplets_j.emplace_back(i+NumUnE,col4_a+2, Vec4_a.operator[](2));//
			}
			triplets_j.emplace_back(i+NumUnE,(3*NumV),((-1)*weight_b[i]*((Dihedral_e1[hedgeid[iJj]])-HEedge[hedgeid[iJj]].dihedral_angle)));
			i++;
			//}
		}
	}
	Jaco_j.setFromTriplets(triplets_j.begin(),triplets_j.end());//Jaco-J

	Jaco_u.resize(3*NumV+1);//3v+k
	for(int j=0;j<NumV;)//,u=0,u<(3*NumV)
	{
		if (!(HEvert_curr[j].isCon))
		{
			
			//cout<<"Jaco_u  "<<"vertIndex "<<j<<endl;
			HE_Edge *out=HEvert[j].outedge;
			HE_Edge *curr=HEvert[j].outedge;
			Vec3d temp=Vec3d(0,0,0);
			do
			{
				Vec3d xi=HEvert_curr[curr->pre->startvert->vertindex].VertPosition;
				Vec3d xj=HEvert_curr[curr->next->startvert->vertindex].VertPosition;////
				temp=temp+cross(xi,xj);
				curr=curr->pair->next;
			}while(curr!=out);
			Jaco_u(3*j)=weight_v*(temp.operator[](0));
			Jaco_u(3*j+1)=weight_v*(temp.operator[](1));
			Jaco_u(3*j+2)=weight_v*(temp.operator[](2));
			//u=u+3;
		}
		else
		{
			Jaco_u(3*j)=0;
			Jaco_u(3*j+1)=0;
			Jaco_u(3*j+2)=0;
		}
		j++;
	}
	Jaco_u(3*NumV)=(-1)*weight_v*(volume_e1-volume);//
}

SparseMatrix<double> A;
//MatrixXd Jaco;
VectorXd b;
VectorXd y;
VectorXd z;
VectorXd delta;

void Gauss_Newton()
{
	A.resize(3*NumV+1,3*NumV+1);
	b.resize(3*NumV+1);
	delta.resize(3*NumV+1);

	double itera_h=1;//#########
	int i=1;
	double E_temp=0;

	VectorXd E;
	VectorXd E_delta;
	E.resize(1);
	E_delta.resize(1);

	HEfaceCurr();
	DihedralCurr();
	Resi_Func();
	Jaco_ju();

	E=0.5*(ResiFunc.transpose()*ResiFunc);
	std::cout<<"E   "<<E<<endl;
	
	//compute delta
	A=Jaco_j.transpose()*Jaco_j;
	b=-1*(Jaco_j.transpose()*ResiFunc.block(0,0,2*NumUnE,1)+Jaco_u*ResiFunc.block(2*NumUnE,0,1,1));
	MatrixXd B=A;
	LDLT<MatrixXd> ldlt;
	ldlt.compute(B);
	y=ldlt.solve(b);
	z=ldlt.solve(Jaco_u);
	delta=y-((z*(Jaco_u.transpose()*y))/(1+(Jaco_u.transpose()*z)));

	do
	{
		for(int ivert=0,iDelta=0;ivert<NumV,iDelta<3*NumV;)
		{
			HEvert_curr[ivert].VertPosition.operator[](0)+=(delta(iDelta)*itera_h);
			HEvert_curr[ivert].VertPosition.operator[](1)+=(delta(iDelta+1)*itera_h);
			HEvert_curr[ivert].VertPosition.operator[](2)+=(delta(iDelta+2)*itera_h);
			iDelta=iDelta+3;
			ivert++;
			//i++;
		} 

		for(int ivert=0;ivert<NumV;ivert++)//initialize the vert_picking
		{
			vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
		}

		interp1+=delta(3*NumV)*itera_h;

		HEfaceCurr();
		DihedralCurr();
		Resi_Func();

		E_delta=0.5*(ResiFunc.transpose()*ResiFunc);
	
		if(E_delta(0)<=E(0))///#######
		{
			break;
		}
		else
		{
			itera_h=0.5*itera_h;//#######
			E(0)=E_delta(0);

			i++;
		}
	} while (itera_h>pow(10,-10));//#######
	std::cout<<i<<endl;
	std::cout<<"final E_delta    "<<E_delta<<endl;
	std::cout<<"final h    "<<itera_h<<endl;
	std::cout<<"final interp1    "<<interp1<<endl;
}


//=====================================================OpenGl========================================================

static double x_offset = 0.0;
static double y_offset = 0.0;
static double x_angle = 0.0;
static double y_angle = 0.0;
static double ScaleSize = 1.0;
static int press_x=0;
static int press_y=0;
static int Mode=0;
int selectVid;
GLuint selectBuffer[128] = {0};
GLint hits;
const GLfloat blue[]= { 0, 0, 1};
const GLfloat red[]= { 1, 0, 0};
int ModelOption=0;
int SelectMode=0, DragStart=0;
vector<int> MultiSelectV;

Vertex MultiCenter;//

#ifndef Box1
#define Box1 1
#endif//
#ifndef Box2
#define Box2 2
#endif//
#ifndef Box3
#define Box3 3
#endif//
#ifndef Cylinder1
#define Cylinder1 4
#endif//
#ifndef Cylinder2
#define Cylinder2 5
#endif//
#ifndef Cylinder3
#define Cylinder3 6//
#endif//
#ifndef Pyramid1
#define Pyramid1 7
#endif//
#ifndef Pyramid2
#define Pyramid2 8
#endif//
#ifndef Pyramid3
#define Pyramid3 9//
#endif//


//***********************************

void InitializeSelect()
{
	SelectMode=0;
	MultiMode=0;
	DragStart=0;
	for (int i = 0; i < MultiSelectV.size(); i++)
	{
		vert_picking[MultiSelectV[i]].selected=false;
	}//
	vert_picking[selectVid].selected=false;
	selectVid=0;
	MultiSelectV.clear();
}


//***********************************

void Draw_Model()
{
	for(int iFace=0;iFace<NumF;iFace++)
	{
		glBegin(GL_TRIANGLES);//GL_TRIANGLES  GL_LINE_LOOP GL_POINTS

			glNormal3d(HEvert_curr[HEface[iFace].firstedge->startvert->vertindex].VertNorm.operator[](0),HEvert_curr[HEface[iFace].firstedge->startvert->vertindex].VertNorm.operator[](1),HEvert_curr[HEface[iFace].firstedge->startvert->vertindex].VertNorm.operator[](2));
			glVertex3d(HEvert_curr[HEface[iFace].firstedge->startvert->vertindex].VertPosition.operator[](0),HEvert_curr[HEface[iFace].firstedge->startvert->vertindex].VertPosition.operator[](1),HEvert_curr[HEface[iFace].firstedge->startvert->vertindex].VertPosition.operator[](2));
			
			glNormal3d(HEvert_curr[HEface[iFace].firstedge->next->startvert->vertindex].VertNorm.operator[](0),HEvert_curr[HEface[iFace].firstedge->next->startvert->vertindex].VertNorm.operator[](1),HEvert_curr[HEface[iFace].firstedge->next->startvert->vertindex].VertNorm.operator[](2));
			glVertex3d(HEvert_curr[HEface[iFace].firstedge->next->startvert->vertindex].VertPosition.operator[](0),HEvert_curr[HEface[iFace].firstedge->next->startvert->vertindex].VertPosition.operator[](1),HEvert_curr[HEface[iFace].firstedge->next->startvert->vertindex].VertPosition.operator[](2));
			
			glNormal3d(HEvert_curr[HEface[iFace].firstedge->pre->startvert->vertindex].VertNorm.operator[](0),HEvert_curr[HEface[iFace].firstedge->pre->startvert->vertindex].VertNorm.operator[](1),HEvert_curr[HEface[iFace].firstedge->pre->startvert->vertindex].VertNorm.operator[](2));
			glVertex3d(HEvert_curr[HEface[iFace].firstedge->pre->startvert->vertindex].VertPosition.operator[](0),HEvert_curr[HEface[iFace].firstedge->pre->startvert->vertindex].VertPosition.operator[](1),HEvert_curr[HEface[iFace].firstedge->pre->startvert->vertindex].VertPosition.operator[](2));
		glEnd();
	}
}

void Draw_Example(vector<HE_Vert> &vert)
{
	for(int iFace=0;iFace<NumF;iFace++)
	{
		glBegin(GL_LINE_LOOP);//GL_TRIANGLES GL_LINE_LOOP GL_POINTS
			glColor3d(0.0,0.0,1.0);

			glNormal3d(vert[HEface[iFace].firstedge->startvert->vertindex].VertNorm.operator[](0),vert[HEface[iFace].firstedge->startvert->vertindex].VertNorm.operator[](1),vert[HEface[iFace].firstedge->startvert->vertindex].VertNorm.operator[](2));
			glVertex3d(vert[HEface[iFace].firstedge->startvert->vertindex].VertPosition.operator[](0),vert[HEface[iFace].firstedge->startvert->vertindex].VertPosition.operator[](1),vert[HEface[iFace].firstedge->startvert->vertindex].VertPosition.operator[](2));

			glNormal3d(vert[HEface[iFace].firstedge->next->startvert->vertindex].VertNorm.operator[](0),vert[HEface[iFace].firstedge->next->startvert->vertindex].VertNorm.operator[](1),vert[HEface[iFace].firstedge->next->startvert->vertindex].VertNorm.operator[](2));
			glVertex3d(vert[HEface[iFace].firstedge->next->startvert->vertindex].VertPosition.operator[](0),vert[HEface[iFace].firstedge->next->startvert->vertindex].VertPosition.operator[](1),vert[HEface[iFace].firstedge->next->startvert->vertindex].VertPosition.operator[](2));
			
			glNormal3d(vert[HEface[iFace].firstedge->pre->startvert->vertindex].VertNorm.operator[](0),vert[HEface[iFace].firstedge->pre->startvert->vertindex].VertNorm.operator[](1),vert[HEface[iFace].firstedge->pre->startvert->vertindex].VertNorm.operator[](2));
			glVertex3d(vert[HEface[iFace].firstedge->pre->startvert->vertindex].VertPosition.operator[](0),vert[HEface[iFace].firstedge->pre->startvert->vertindex].VertPosition.operator[](1),vert[HEface[iFace].firstedge->pre->startvert->vertindex].VertPosition.operator[](2));
		glEnd();
	}
}

void Draw_Vert(GLenum mode)
{
	glPushMatrix();
	for(int ivert=0;ivert<NumV;ivert++)
	{
		if(mode==GL_SELECT)
		{
			glLoadName(ivert);
		}
		vert_picking[ivert].Draw();
	}
	glPopMatrix();
}


//***********************************

void ProcessMenuEvents(int option) 
{
	 switch (option) 
	 {
		 case Box1 :
			 {
				cout<<"box and example1"<<endl;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				ModelOption=Box1;
				MESH= new ObjMesh("1 box-same.obj");
				example_1= new ObjMesh("1-box-same-stretch x-y-.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(4);
				HEvert_Con[0]=0;
				HEvert_Con[1]=1;
				HEvert_Con[2]=6;
				HEvert_Con[3]=7;//不动的点
				HalfEdge_MESH();
				HalfEdge_example_1();
				volume_e1=Volume(example_1,HEvert_e1);//
			
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector不动点的边
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
         case Box2 : 
			 {
				cout<<"box and example2"<<endl; 
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				ModelOption=Box2;
				MESH= new ObjMesh("box.obj");
				example_1= new ObjMesh("box_e2.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(4);
				HEvert_Con[0]=0;
				HEvert_Con[1]=1;
				HEvert_Con[2]=6;
				HEvert_Con[3]=7;//
				HalfEdge_MESH();
				HalfEdge_example_1();
				volume_e1=Volume(example_1,HEvert_e1);//
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
		case Box3 : 
			{
				cout<<"box and two examples"<<endl;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				ModelOption=Box3;
				MESH= new ObjMesh("box.obj");
				example_1= new ObjMesh("box_e1.obj");
				example_2= new ObjMesh("box_e2.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(4);
				HEvert_Con[0]=0;
				HEvert_Con[1]=1;
				HEvert_Con[2]=6;
				HEvert_Con[3]=7;//
				HalfEdge_MESH();
				HalfEdge_example_1();
				HalfEdge_example_2();
				volume_e1=Volume(example_1,HEvert_e1);
				volume_e2=Volume(example_2,HEvert_e2);//
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			}
         case Cylinder1 : 
			 {
 				cout<<"cylinder and example1"<<endl; 
				ModelOption=Cylinder1;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				MESH= new ObjMesh("2.1 cylinder box.obj");
				example_1= new ObjMesh("2.1 cylinder box-z-90.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(11);
				for (int ivert = 0; ivert < 4; ivert++)
				{
					HEvert_Con[ivert]=ivert;
				}
				//HEvert_Con[10]=20;//
				HalfEdge_MESH();
				HalfEdge_example_1();
				volume_e1=Volume(example_1,HEvert_e1);//
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();		
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
         case Cylinder2 : 
			 {
				cout<<"cylinder and example2"<<endl;
				ModelOption=Cylinder2;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				MESH= new ObjMesh("7_s.obj");
				example_1= new ObjMesh("60_s.obj");//
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(4);
				for (int ivert = 0; ivert < 4; ivert++)
				{
					HEvert_Con[ivert]=ivert;
				}//
				HalfEdge_MESH();
				HalfEdge_example_1();
				volume_e1=Volume(example_1,HEvert_e1);//
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
		 case Cylinder3 : 
			 {
				cout<<"cylinder and two examples"<<endl; 
				ModelOption=Cylinder3;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				MESH= new ObjMesh("cylinder.obj");
				example_1= new ObjMesh("cylinder_e1.obj");
				example_2= new ObjMesh("cylinder_e2.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(11);
				for (int ivert = 0; ivert < 10; ivert++)
				{
					HEvert_Con[ivert]=ivert;
				}//
				HEvert_Con[10]=20;
				HalfEdge_MESH();
				HalfEdge_example_1();
				HalfEdge_example_2();
				volume_e1=Volume(example_1,HEvert_e1);
				volume_e2=Volume(example_2,HEvert_e2);//
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
		 case Pyramid1 : 
			 {
				cout<<"pyramid and example1"<<endl;  
				ModelOption=Pyramid1;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				MESH= new ObjMesh("pyramid0.obj");
				example_1= new ObjMesh("pyramid0_e1.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(2);
				HEvert_Con[0]=0;
				HEvert_Con[1]=2;
				HalfEdge_MESH();
				HalfEdge_example_1();
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
         case Pyramid2 : 
			 {
				cout<<"pyramid and example2"<<endl;
				ModelOption=Pyramid2;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				MESH= new ObjMesh("pyramid0.obj");
				example_1= new ObjMesh("pyramid0_e2.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(2);
				HEvert_Con[0]=0;
				HEvert_Con[1]=2;
				HalfEdge_MESH();
				HalfEdge_example_1();
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
			 
		 case Pyramid3 : 
			 {
				cout<<"pyramid and two examples"<<endl; 
				ModelOption=Pyramid3;
				HEvert.clear();
				HEface.clear();
				HEedge.clear();//
				HEvert_curr.clear();
				HEface_curr.clear();
				Dihedral_curr.clear();//
				HEvert_e1.clear();
				HEface_e1.clear();
				Dihedral_e1.clear();//
				InitializeSelect();//
				MESH= new ObjMesh("pyramid0.obj");
				example_1= new ObjMesh("pyramid0_e1.obj");
				example_2= new ObjMesh("pyramid0_e2.obj");
				NumV = MESH->getNumVertices();
				NumF = MESH->getNumFaces();
				NumG = MESH->getNumGroups(); 
				Edgelist.clear();
				Edgelist.resize(NumV);
				for(int iCol=0;iCol<NumV;iCol++)
				{
					Edgelist[iCol].resize(NumV);
				}//
				HEvert_Con.clear();
				HEvert_Con.resize(2);
				HEvert_Con[0]=0;
				HEvert_Con[1]=2;
				HalfEdge_MESH();
				HalfEdge_example_1();
				HalfEdge_example_2();
				HEvertCurr();
				HEfaceCurr();
				DihedralCurr();
				HEvnormCurr();
				for (int ivert = 0; ivert < HEvert_Con.size(); ivert++)
				{
					for (int jvert = ivert+1; jvert < HEvert_Con.size(); jvert++)
					{
						if(Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].isEdge==true)
						{Edgelist[HEvert_Con[ivert]][HEvert_Con[jvert]].HEedgeID=-1;}
					}
				}
				hedgeid.clear();
				hedgeid.resize(1.5*NumF);
				int iedge=0;
				for (int irow= 0; irow < NumV; irow++)//Edge index
				{
					for(int icol=irow;icol<NumV;icol++)
					{
						if (Edgelist[irow][icol].isEdge)
						{
							hedgeid[iedge]=Edgelist[irow][icol].HEedgeID;
							iedge++;
						}
					}
				}// two for-loop for hedgeid vector
				vert_picking.clear();
				vert_picking.resize(NumV);
				for(int ivert=0;ivert<NumV;ivert++)
				{
					vert_picking[ivert].vertpos=HEvert_curr[ivert].VertPosition;
					vert_picking[ivert].iscon=HEvert_curr[ivert].isCon;
				}//for-loop to initialize the vert_picking
				weight_s.clear();
				weight_b.clear();
				Weight();
				NumUnE=weight_s.size();//
				triplets_j.clear();
				glutPostRedisplay();
				break;
			 }
	 }
}

void CreateGLUTMenus() 
{
         int menu,submenu_box,submenu_cylinder,submenu_pyramid;
 
         submenu_box = glutCreateMenu(ProcessMenuEvents);
		 glutAddMenuEntry("Box + example1 ",Box1);
		 glutAddMenuEntry("Box + example2 ",Box2);
		 glutAddMenuEntry("Box + example1 + example2 ",Box3);

		 submenu_cylinder = glutCreateMenu(ProcessMenuEvents);
		 glutAddMenuEntry("Cylinder + example1",Cylinder1);
		 glutAddMenuEntry("Cylinder + example2",Cylinder2);
		 glutAddMenuEntry("Cylinder + example1 + example2 ",Cylinder3);

		 submenu_pyramid = glutCreateMenu(ProcessMenuEvents);
		 glutAddMenuEntry("Pyramid + example1",Pyramid1);
		 glutAddMenuEntry("Pyramid + example2",Pyramid2);
		 glutAddMenuEntry("Pyramid + example1 + example2 ",Pyramid3);
 
         menu = glutCreateMenu(ProcessMenuEvents);
		 glutAddSubMenu("Box",submenu_box);
		 glutAddSubMenu("Cylinder",submenu_cylinder);
		 glutAddSubMenu("Pyramid",submenu_pyramid);
         glutAttachMenu(GLUT_RIGHT_BUTTON);
}


//***********************************

void PressSpecialKeys(int key, int x, int y) 
{    
	if (key==GLUT_KEY_UP)
	{
		double OldSize = ScaleSize;
		ScaleSize =ScaleSize*1.5;
		glutPostRedisplay();
	}
	if (key==GLUT_KEY_DOWN)
	{
		double OldSize = ScaleSize;
		ScaleSize =ScaleSize*0.6;
		glutPostRedisplay();
	}
} 


//***********************************

void SelectObjects(int x,int y)//x,y是当前鼠标的坐标
{
	GLint viewport[4];//存放可视区参数
		
	glGetIntegerv(GL_VIEWPORT,viewport);//获取当前视口坐标参数
	glSelectBuffer(128,selectBuffer);//选择名称栈存放被选中的名称

	glRenderMode(GL_SELECT);//设置当前为 选择模式
	//初始化名称栈
	glInitNames();
	glPushName(-1);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
		glLoadIdentity();
		gluPickMatrix(x,viewport[3]-y,10.0, 10.0, viewport);//创建用于选择的投影矩阵栈
		gluPerspective(60, 1, 0.1, 100);
		glMatrixMode(GL_MODELVIEW);
		Draw_Vert(GL_SELECT);
	glPopMatrix();
	
	glFlush();
	
	hits = glRenderMode(GL_RENDER);
	
	glMatrixMode(GL_PROJECTION);  
	glLoadIdentity();  
	glViewport(0, 0, 1000, 1000);  
	gluPerspective(60, 1, 0.1, 100);

	if (hits>0)
	{
		int n=0;
		double minz=selectBuffer[1];//0x7fffffff
		for(int i=1;i<hits;i++)
		{
			if ((selectBuffer[1+i*4])<minz)//0x7fffffff
			{
				n=i;
				minz=selectBuffer[1+i*4];//0x7fffffff;
			}
		}
		vert_picking[selectBuffer[3+n*4]].selected=!(vert_picking[selectBuffer[3+n*4]].selected);//;true
		selectVid=selectBuffer[3+n*4];
	}
}

void Drag(int x,int y)
{
	if(hits>0)
	{
		GLint viewport[4];//存放可视区参数
		GLdouble  modelview[16], projection[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetIntegerv(GL_VIEWPORT,viewport);

		GLdouble wVx,wVy,wVz;
		GLdouble Vx,Vy,Vz;
		double y_new=viewport[3]-y;

		if (MultiMode>0)
		{
			Vec3d MultiCenter_temp=MultiCenter.vertpos;
			gluProject(MultiCenter.vertpos.operator[](0),MultiCenter.vertpos.operator[](1),MultiCenter.vertpos.operator[](2),
				       modelview,projection,viewport,&wVx, &wVy, &wVz);//selected vert 3D coord change to 2D screen coord
			gluUnProject(x,y_new,wVz,modelview,projection,viewport,&Vx, &Vy, &Vz);//mouse screen coord changes to 3D coord

			//MultiCenter.vertpos.operator[](0)+=(0.8); 
			//MultiCenter.vertpos.operator[](1)+=(0.4);
			//MultiCenter.vertpos.operator[](2)+=0;//for test

			Vec3d replacement=MultiCenter.vertpos-MultiCenter_temp;
			cout<<"replacement   "<<replacement<<endl;
			for (int i = 0; i < MultiSelectV.size(); i++)
			{
				vert_picking[MultiSelectV[i]].vertpos=vert_picking[MultiSelectV[i]].vertpos+replacement;//(replacement/MultiSelectV.size());
				HEvert_curr[MultiSelectV[i]].VertPosition=vert_picking[MultiSelectV[i]].vertpos;
				cout<<"MultiSelectV:     "<<MultiSelectV[i]<<endl;
			}
		}
		else
		{
			cout<<"before dragging HEvert_curr:"<<"\t";
			cout<<selectVid<<"\t"<<HEvert_curr[selectVid].VertPosition<<endl;
			gluProject(vert_picking[selectVid].vertpos.operator[](0), 
				   vert_picking[selectVid].vertpos.operator[](1), 
				   vert_picking[selectVid].vertpos.operator[](2),
				   modelview,projection,viewport,&wVx, &wVy, &wVz);//selected vert 3D coord change to 2D screen coord
			gluUnProject(x,y_new,wVz,modelview,projection,viewport,&Vx, &Vy, &Vz);//mouse screen coord changes to 3D coord

			//vert_picking[selectVid].vertpos.operator[](0)+=(0); 
			//vert_picking[selectVid].vertpos.operator[](1)+=0;
			//vert_picking[selectVid].vertpos.operator[](2)+=0;//for test

			HEvert_curr[selectVid].VertPosition=vert_picking[selectVid].vertpos;

			cout<<"after dragging HEvert_curr:"<<"\t";
			cout<<selectVid<<"\t"<<HEvert_curr[selectVid].VertPosition<<endl;
			cout<<"after dragging vert_picking:"<<"\t";
			cout<<selectVid<<"\t"<<vert_picking[selectVid].vertpos<<endl;
		}

		Gauss_Newton();
	}
	else
	{
		cout<<"No picking"<<endl;
	}
	
}

void PressKeys(unsigned char key, int x, int y) 
{    
	if (key==27)
	{
		exit(0);
	}
	if (key=='s'||key=='S')//
	{
		SelectMode=1;
		MultiMode=1;
	}
	//=glutGetModifiers();
} 

void ReleaseKeys(unsigned char key, int x, int y) 
{    
	if (key=='s'||key=='S')//
	{
		SelectMode=0;
		DragStart=1;

		Vec3d verpos_temp=Vec3d(0,0,0);
		sort(MultiSelectV.begin(),MultiSelectV.end());
		MultiSelectV.erase( unique( MultiSelectV.begin(), MultiSelectV.end() ), MultiSelectV.end() );//delete same element
		for (int i = 0; i < MultiSelectV.size(); i++)
		{
			verpos_temp+=vert_picking[MultiSelectV[i]].vertpos;
		}
		MultiCenter.vertpos=verpos_temp/MultiSelectV.size();
		cout<<"MultiCenter   "<<MultiCenter.vertpos<<endl;
		MultiCenter.selected=true;//!MultiCenter.selected
		//cout<<"DragStart:  "<<DragStart<<endl;
	}
	glutPostRedisplay();
} 


//***********************************

void Display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  
	glEnable(GL_DEPTH_TEST);  
	glEnable(GL_LIGHTING);  
	glEnable(GL_LIGHT0);  

	glClearColor(0,0,0,1.0);  
	glClearDepth(1.0);  
	glDepthFunc(GL_LESS);  
	glEnable(GL_DEPTH_TEST);
    
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, 0.1, 100);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,1,6,0,1,0,0,1,0);

	glRotated(x_angle, 0.0, 1.0, 0.0);
	glRotated(y_angle, 1.0, 0.0, 0.0);
	glScaled(ScaleSize, ScaleSize, ScaleSize);
	
	//=============================================
	if (ModelOption==0)
	{
		MESH= new ObjMesh("box.obj");
		NumV = MESH->getNumVertices();
		NumF = MESH->getNumFaces();
		NumG = MESH->getNumGroups(); 
		Edgelist.resize(NumV);
		for(int iCol=0;iCol<NumV;iCol++)
		{
			Edgelist[iCol].resize(NumV);
		}//
		HalfEdge_MESH();

		vert_picking.resize(NumV);
		HEvert_curr.resize(NumV);
		for(int ivert=0;ivert<NumV;ivert++)
		{
			vert_picking[ivert].vertpos=HEvert[ivert].VertPosition;
			HEvert_curr[ivert].VertPosition=HEvert[ivert].VertPosition;
		}//for-loop to initialize the vert_picking
	}//#######
	else if (ModelOption==Box1)
	{
		glPushMatrix();
			glTranslated(0,2,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();
	}
	else if (ModelOption==Box2)
	{
		glPushMatrix();
			glTranslated(0,2,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();
	}
	else if (ModelOption==Box3)
	{
		glPushMatrix();
			glTranslated(-1,2,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();//
		glPushMatrix();
			glTranslated(1,2,0);
			Draw_Example(HEvert_e2);
		glPopMatrix();//
	}/////////////Box//////////////
	else if (ModelOption==Cylinder1)
	{
		Draw_Example(HEvert_e1);
	}
	else if (ModelOption==Cylinder2)
	{
		glPushMatrix();
			glTranslated(0,2.5,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();//
	}
	else if (ModelOption==Cylinder3)
	{
		glPushMatrix();
			glTranslated(-1,2.5,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();//
		glPushMatrix();
			glTranslated(1,2.5,0);
			Draw_Example(HEvert_e2);
		glPopMatrix();//
	}/////////Cylinder/////////
	else if (ModelOption==Pyramid1)
	{
		glPushMatrix();
			glTranslated(0,2,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();//
	}
	else if (ModelOption==Pyramid2)
	{
		glPushMatrix();
			glTranslated(0,2,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();//
	}
	else if (ModelOption==Pyramid3)
	{
		glPushMatrix();
			glTranslated(-1,2,0);
			Draw_Example(HEvert_e1);
		glPopMatrix();//
		glPushMatrix();
			glTranslated(1,2,0);
			Draw_Example(HEvert_e2);
		glPopMatrix();//
	}/////////Pyramid/////////


	Draw_Vert(GL_RENDER);//Loadname

	if (MultiMode>0)
	{
		MultiCenter.Draw();

	}

	if (ModelOption==0)
	{
		glColor4d(1.0,1.0,1.0,0.4);
	}
	else
	{
		glColor4d(1.0,1.0,0.0,0.4);
	}
	Draw_Model();

	glDisable(GL_LIGHT0);  
	glDisable(GL_LIGHTING);  
	glDisable(GL_DEPTH_TEST); 
	glutSwapBuffers();
}

void Reshape(int w, int h)
{
	glViewport(0, 0, 1000, 1000);  // 注意默认的reshape函数中的glViewport函数设置 
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	gluPerspective(60, 1, 0.1, 100);
}

void Mouse(int button, int state, int mousex, int mousey) 
{
	if(state == GLUT_UP)
	{ 
		return;
	}  
	else
	{
		if (state == GLUT_DOWN)
		{
			press_x=mousex;
			press_y=mousey;
			if (button == GLUT_LEFT_BUTTON)
			{
				Mode =1;//select points
				if((SelectMode>0)&&(MultiMode>0))//(MultiMode>0)
				{
					SelectObjects(mousex,mousey);
					MultiSelectV.emplace_back(selectVid);
				}
				else if ((MultiMode>0)&&(DragStart>0))
				{
					Drag(mousex,mousey);
					cout<<"HEvert after deformation:  "<<endl;
					for (int ivert = 0; ivert < NumV; ivert++)
					{
						cout<<ivert<<"\t"<<HEvert[ivert].VertPosition<<"\t\t"<<HEvert_curr[ivert].VertPosition<<endl;
					}
					cout<<"############################"<<endl;
				}

				else if((SelectMode==0)&&(MultiMode==0))
				{
					vert_picking[selectVid].selected=false;
					SelectObjects(mousex,mousey);
					Drag(mousex,mousey);
					cout<<"HEvert after deformation:  "<<endl;
					for (int ivert = 0; ivert < NumV; ivert++)
					{
						cout<<ivert<<"\t"<<HEvert[ivert].VertPosition<<"\t\t"<<HEvert_curr[ivert].VertPosition<<endl;
					}
					cout<<"############################"<<endl;
				}
			}
			else if (button == GLUT_MIDDLE_BUTTON)
			{
				Mode = 2;
				cout<<"Rotate"<<endl;//Rotate
			}
			else if (button == GLUT_RIGHT_BUTTON)
			{
				Mode = 3;
				cout<<"Scale"<<endl;//Scale
			}
		}
	}
		glutPostRedisplay();
}

void Motion(int x, int y)
	{
		if (Mode == 1)
		{
			press_x=x;
			press_y=y;
		}////
		else if(Mode==2)
		{
			x_angle += (x - press_x) / 5.0;
			if (x_angle > 180)
			{
				x_angle -= 360;
			}
			else if (x_angle <-180)
			{
				x_angle += 360;
			}
			press_x = x;//

			y_angle += (y - press_y) / 5.0;
			if (y_angle > 180)
			{
				y_angle -= 360;
			}
			else if (y_angle <-180)
			{
				y_angle += 360;
			}
			press_y = y;//
		}////
		else if (Mode == 3)
		{
			double OldSize = ScaleSize;
			ScaleSize =ScaleSize*(1+((y - press_y)/100.0));
			press_y=y;
		}////
		glutPostRedisplay(); 
	}//mouse motion

void Lights()  
{  
    GLfloat AmbientLight[]  ={0.5f,  0.5f,  0.5f,  1.0f};
    GLfloat DiffuseLight[]  ={0.9f,  0.9f,  0.9f,  1.0f};  
    GLfloat SpecularLight[] ={1.0f,  1.0f,  1.0f,  1.0f};
    GLfloat LightPos[]={1.0f,1.0f,5.0f,0.0f};
	
	glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);            
    glLightfv(GL_LIGHT0, GL_AMBIENT,  AmbientLight);    
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  DiffuseLight);      
    glLightfv(GL_LIGHT0, GL_SPECULAR, SpecularLight);   
    glLightfv(GL_LIGHT0, GL_POSITION, LightPos);          
    glEnable(GL_LIGHT0);            

	glEnable(GL_COLOR_MATERIAL);    
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);  
    glMaterialfv(GL_FRONT, GL_SPECULAR, SpecularLight); 
    glMateriali(GL_FRONT, GL_SHININESS, 100);         
}


//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

int main(int argc, char **argv)
{
	//===========================================OpenGL======================================================== 
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(1000,1000);
	glutInitWindowPosition(100,100);
	glutCreateWindow("Mesh Deformation");
	
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);
	
	Lights();
	
	CreateGLUTMenus();
	
	glutKeyboardFunc(PressKeys);
	glutKeyboardUpFunc(ReleaseKeys);
	glutSpecialFunc(PressSpecialKeys);

	glutMainLoop();

	return 1;
}