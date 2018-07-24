#include "VECTOR.h"
#include "Face.h"
#include "Mesh.h"
#include "Init.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>

using namespace std;
float cameraAngle = 0.0f;
float cameraPos = 90.0f;
float cameraX = 0.0f, cameraY = 30.0f, cameraZ = 90.0f;

VECTOR3D eye = VECTOR3D(0.0f, 30.0f, 90.0f);
VECTOR3D light = { 0.0, 110.0, 0.0 };

int depth = 3;
int isoctree = 0;
int start = 0;
int treesize = 50;

clock_t startTime, finishTime;

GLdouble mvMatrix[16];
GLdouble projMatrix[16];
GLint viewport[4];

struct Ray
{
	//   ray의 원점
	VECTOR3D origin;
	//   ray의 방향
	VECTOR3D dir;

	Ray(VECTOR3D origin, VECTOR3D dir)
	{
		this->origin = origin;
		this->dir = dir;
	}
};

class object
{
public:
	int mode;
	VECTOR3D k_ambient;
	VECTOR3D k_diffuse;
	VECTOR3D k_specular;
	float k_shineness;

	object() {}

	virtual ~object() {}

	virtual bool hit(Ray r, float *t) = 0 {}
	virtual VECTOR3D getColor(VECTOR3D point, VECTOR3D light, VECTOR3D ray) = 0 {}
	virtual VECTOR3D get_normal(VECTOR3D point) = 0 {};
	virtual int get_mode() = 0 {};
	virtual bool exist(VECTOR3D min, VECTOR3D max) = 0 {}
};

class sphere : public object
{
public:
	int mode;
	VECTOR3D cen;
	float rad;

	VECTOR3D k_ambient;
	VECTOR3D k_diffuse;
	VECTOR3D k_specular;
	float k_shineness;

	virtual bool hit(Ray r, float *t) {
		float b = 2.0 * (r.dir.InnerProduct(r.origin - cen));
		float c = pow((r.origin - cen).Magnitude(), 2) - (rad*rad);

		float d = b * b - (4.0*c);   //판별식

		if (d < 0.0)
		{
			*t = 10000.0;
			return false;
		}

		else if (d == 0.0)
		{
			*t = (-1.0*b) / 2.0;

			if (*t > 0.01)
				return true;
			else
				return false;
		}

		else
		{
			float t1 = ((-1.0*b) - sqrt(d)) / 2.0;
			float t2 = ((-1.0 *b) + sqrt(d)) / 2.0;

			if (t1 > 0.01 && t2 > 0.01)
			{
				*t = (t1 < t2) ? t1 : t2;
				return true;
			}
			else
				return false;
		}
	}

	virtual VECTOR3D getColor(VECTOR3D point, VECTOR3D light, VECTOR3D ray_origin)
	{
		VECTOR3D N = get_normal(point);
		N.Normalize();

		VECTOR3D I = ray_origin - point;
		I.Normalize();

		VECTOR3D L = light - point;
		L.Normalize();

		VECTOR3D R = (-1.0 * L) + 2.0 * ((L.InnerProduct(N)) * N);

		float diffuse = std::max((float)0.0, N.InnerProduct(L));
		float specular = pow(std::max((float)0.0, I.InnerProduct(R)), k_shineness);

		return diffuse * k_diffuse + specular * k_specular;
	}

	virtual VECTOR3D get_normal(VECTOR3D point)
	{
		return point - cen;
	}

	virtual int get_mode()
	{
		return mode;
	}

	virtual bool exist(VECTOR3D min, VECTOR3D max)
	{
		VECTOR3D o = { (max.x + min.x) / 2, (max.y + min.y) / 2, (max.z + min.z) / 2 };

		float oc = (o - cen).Magnitude();
		float mm = (max - min).Magnitude() / 2 + rad;
		if (oc <= mm)
			return true;
		else
			return false;
	}
};

class triangle : public object
{
public:
	int mode;
	VECTOR3D a;
	VECTOR3D b;
	VECTOR3D c;

	VECTOR3D va;
	VECTOR3D vb;
	VECTOR3D vc;

	VECTOR3D k_ambient;
	VECTOR3D k_diffuse;
	VECTOR3D k_specular;
	float k_shineness;


	virtual bool hit(Ray r, float *t)
	{
		VECTOR3D v_ba = b - a;
		VECTOR3D v_ca = c - a;
		VECTOR3D pvec = r.dir.CrossProduct(v_ca);
		float det = v_ba.InnerProduct(pvec);

		float invDet = 1 / det;
		VECTOR3D tvec = r.origin - a;
		float beta = tvec.InnerProduct(pvec) * invDet;
		VECTOR3D qvec = tvec.CrossProduct(v_ba);
		float gamma = r.dir.InnerProduct(qvec) * invDet;
		float alpha = 1 - beta - gamma;

		VECTOR3D N = (b - a).CrossProduct(c - a);


		if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1 || N.InnerProduct(r.dir) >= 0)
			return false;
		else
		{
			//cout << "tri" << endl;
			float D = N.InnerProduct(a);
			int dis = (-N.InnerProduct(r.origin) + D) / N.InnerProduct(r.dir);

			if (dis <= 0)
				return false;

			*t = dis;
			return true;
		}
	}

	virtual VECTOR3D getColor(VECTOR3D point, VECTOR3D light, VECTOR3D ray_origin)
	{
		VECTOR3D N = get_normal(point);
		N.Normalize();

		VECTOR3D I = ray_origin - point;
		I.Normalize();

		VECTOR3D L = light - point;
		L.Normalize();

		VECTOR3D R = (-1.0 * L) + 2.0 * ((L.InnerProduct(N)) * N);

		float diffuse = std::max((float)0.0, N.InnerProduct(L));
		float specular = pow(std::max((float)0.0, I.InnerProduct(R)), k_shineness);

		return diffuse * k_diffuse + specular * k_specular;
	}

	virtual VECTOR3D get_normal(VECTOR3D point)
	{
		//////////////////////////////////////triangle
		float da = sqrt(pow((point.x - a.x), 2.0) + pow((point.y - a.y), 2.0) + pow((point.z - a.z), 2.0));
		float db = sqrt(pow((point.x - b.x), 2.0) + pow((point.y - b.y), 2.0) + pow((point.z - b.z), 2.0));
		float dc = sqrt(pow((point.x - c.x), 2.0) + pow((point.y - c.y), 2.0) + pow((point.z - c.z), 2.0));

		float dd = da + db + dc;

		float al = (1 - da / dd);
		float be = (1 - db / dd);
		float ca = (1 - dc / dd);

		VECTOR3D pn = (al*va + be * vb + ca * vc);

		return pn;
		
		
		///////////////////////////////////// plane
		//VECTOR3D N = (b - a).CrossProduct(c - a);
		//return N;
	}

	virtual int get_mode()
	{
		return mode;
	}

	virtual bool exist(VECTOR3D min, VECTOR3D max)
	{
		/*

		if (a.x > min.x && a.x < max.x && a.y > min.y && a.y < max.y && a.z > min.z && a.z < max.z
		|| b.x > min.x && b.x < max.x && b.y > min.y && b.y < max.y && b.z > min.z && b.z < max.z
		|| c.x > min.x && c.x < max.x && c.y > min.y && c.y < max.y && c.z > min.z && c.z < max.z)
		{
		return true;
		}
		else
		{*/
		//   return false;
		VECTOR3D N = (b - a).CrossProduct(c - a);
		N.Normalize();

		VECTOR3D o = { (max.x + min.x) / 2, (max.y + min.y) / 2, (max.z + min.z) / 2 };

		float distance = abs((a - o).InnerProduct(N));

		if (distance <= (max - min).Magnitude() / 2)
			return true;
		else
			return false;
		//}
	}
};

vector<object*> objects;

Mesh mesh1 = Mesh();
Mesh mesh2 = Mesh();
Mesh mesh3 = Mesh();
Mesh mesh4 = Mesh();
Mesh mesh5 = Mesh();

void Meshload()
{

	mesh3.LoadMesh("charb.obj");
}

void ComputeNormal()
{
	mesh3.ComputeFaceNormal();
	mesh3.FindNeighborFaceArray();
	mesh3.ComputeVertexNormal();
}


void RenderMesh(Mesh mesh)
{
	for (int i = 0; i < mesh.faceArray.size(); i++)
	{
		Vertex v0 = mesh.vertexArray[mesh.faceArray[i].vertex0];
		Vertex v1 = mesh.vertexArray[mesh.faceArray[i].vertex1];
		Vertex v2 = mesh.vertexArray[mesh.faceArray[i].vertex2];


		//cout << v0.position.x << endl;
		triangle *tmp = new triangle();

		//cout<<&tmp<<endl;

		tmp->mode = 0;
		tmp->a = { v0.position.x, v0.position.y, v0.position.z };

		tmp->va = { v0.normal.x,v0.normal.y,v0.normal.z };
		tmp->b = { v1.position.x, v1.position.y, v1.position.z };
		tmp->vb = { v1.normal.x,v1.normal.y,v1.normal.z };
		tmp->c = { v2.position.x, v2.position.y, v2.position.z };
		tmp->vc = { v2.normal.x,v2.normal.y,v2.normal.z };
		tmp->k_ambient = { 0.3f,0.1f,0.1f };
		tmp->k_diffuse = { 4.0f,1.0f,1.5f };
		tmp->k_specular = { 15.0f,5.0f,10.0f };
		tmp->k_shineness = 25.0f;

		objects.push_back(tmp);
	}
}


#pragma region Octree
struct TreeNode
{
	TreeNode* parent;
	TreeNode* children[8];

	VECTOR3D boxMin;
	VECTOR3D boxMax;

	vector<object*> *objs;

	TreeNode()
	{
		parent = nullptr;
		objs = nullptr;
		for (int i = 0; i < 8; i++)
		{
			children[i] = nullptr;
		}
	}
};

void Subdivide(TreeNode* node, int depth)
{
	if (depth >= 2)
		return;

	for (int i = 0; i < 8; i++)
	{
		node->children[i] = new TreeNode;
		node->children[i]->parent = node;
	}

	//   left top back
	node->children[0]->boxMin.x = node->boxMin.x;
	node->children[0]->boxMax.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[0]->boxMin.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[0]->boxMax.y = node->boxMax.y;
	node->children[0]->boxMin.z = node->boxMin.z;
	node->children[0]->boxMax.z = (node->boxMin.z + node->boxMax.z) / 2;
	//   right top back
	node->children[1]->boxMin.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[1]->boxMax.x = node->boxMax.x;
	node->children[1]->boxMin.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[1]->boxMax.y = node->boxMax.y;
	node->children[1]->boxMin.z = node->boxMin.z;
	node->children[1]->boxMax.z = (node->boxMin.z + node->boxMax.z) / 2;
	//   left bottom back
	node->children[2]->boxMin.x = node->boxMin.x;
	node->children[2]->boxMax.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[2]->boxMin.y = node->boxMin.y;
	node->children[2]->boxMax.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[2]->boxMin.z = node->boxMin.z;
	node->children[2]->boxMax.z = (node->boxMin.z + node->boxMax.z) / 2;
	//   right bottom back
	node->children[3]->boxMin.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[3]->boxMax.x = node->boxMax.x;
	node->children[3]->boxMin.y = node->boxMin.y;
	node->children[3]->boxMax.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[3]->boxMin.z = node->boxMin.z;
	node->children[3]->boxMax.z = (node->boxMin.z + node->boxMax.z) / 2;

	//   left top front
	node->children[4]->boxMin.x = node->boxMin.x;
	node->children[4]->boxMax.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[4]->boxMin.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[4]->boxMax.y = node->boxMax.y;
	node->children[4]->boxMin.z = (node->boxMin.z + node->boxMax.z) / 2;
	node->children[4]->boxMax.z = node->boxMax.z;
	//   right top front
	node->children[5]->boxMin.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[5]->boxMax.x = node->boxMax.x;
	node->children[5]->boxMin.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[5]->boxMax.y = node->boxMax.y;
	node->children[5]->boxMin.z = (node->boxMin.z + node->boxMax.z) / 2;
	node->children[5]->boxMax.z = node->boxMax.z;
	//   left bottom front
	node->children[6]->boxMin.x = node->boxMin.x;
	node->children[6]->boxMax.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[6]->boxMin.y = node->boxMin.y;
	node->children[6]->boxMax.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[6]->boxMin.z = (node->boxMin.z + node->boxMax.z) / 2;
	node->children[6]->boxMax.z = node->boxMax.z;
	//   right bottom front
	node->children[7]->boxMin.x = (node->boxMin.x + node->boxMax.x) / 2;
	node->children[7]->boxMax.x = node->boxMax.x;
	node->children[7]->boxMin.y = node->boxMin.y;
	node->children[7]->boxMax.y = (node->boxMin.y + node->boxMax.y) / 2;
	node->children[7]->boxMin.z = (node->boxMin.z + node->boxMax.z) / 2;
	node->children[7]->boxMax.z = node->boxMax.z;


	for (int i = 0; i < 8; i++)
	{
		Subdivide(node->children[i], depth + 1);
	}
}

void PostOrderTraversal(TreeNode* node)
{
	for (int i = 0; i < 8; i++)
	{
		//   자식 node가 있다면 자식 node 먼저 visit
		if (node->children[i] != nullptr)
			PostOrderTraversal(node->children[i]);

		if (i == 7)
		{
			//   visit
			std::cout << "bounding box min x: " << node->boxMin.x <<
				" y: " << node->boxMin.y << " z: " << node->boxMin.z << std::endl;
			std::cout << "bounding box max x: " << node->boxMax.x <<
				" y: " << node->boxMax.y << " z: " << node->boxMax.z << std::endl;

			if (objects[0]->exist(node->boxMin, node->boxMax))
				cout << "exist" << endl;

			if (node->objs != nullptr)
				std::cout << node->objs->size() << endl;

			std::cout << std::endl;
		}
	}
}
#pragma endregion


//   한 공간에 하나의 오브젝트만 존재할 때까지 공간 분할
void SpaceDivision(TreeNode *node, vector<object*> curObjects)
{
	vector<object*> *newObjs = new vector<object*>();

	for (int i = 0; i < curObjects.size(); i++)
	{
		object *obj = curObjects[i];

		if (obj->exist(node->boxMin, node->boxMax))
		{
			newObjs->push_back(obj);
			//cout << newObjects.size() << endl;
		}
	}

	if (newObjs->size() > treesize)
	{
		Subdivide(node, 1);

		for (int i = 0; i < 8; i++)
			SpaceDivision(node->children[i], *newObjs);
	}
	else if (newObjs->size() > 0)
	{
		//vector<object*> *objs = new vector<object*>();

		//for (int i = 0; i < newObjs.size(); i++)
		//{
		//	//cout << "hello" << endl;
		//	objs->push_back(newObjs[i]);
		//}
		//cout << node->objs->size() << endl;

		node->objs = newObjs;
	}
}


//   ray가 공간과 교차하는지 검사
bool RayTraversal(TreeNode* octree, Ray ray1) {
	Ray ray = ray1;
	if (ray.dir.x == 0)
		return true;

	// fixes for rays with negative direction
	if (ray.dir.x < 0) {
		ray.origin.x = octree->boxMin.x + octree->boxMax.x - ray.origin.x;
		ray.dir.x = -ray.dir.x;
	}
	if (ray.dir.y < 0) {
		ray.origin.y = octree->boxMin.y + octree->boxMax.y - ray.origin.y;
		ray.dir.y = -ray.dir.y;
	}
	if (ray.dir.z < 0) {
		ray.origin.z = octree->boxMin.z + octree->boxMax.z - ray.origin.z;
		ray.dir.z = -ray.dir.z;
	}

	double divx = 1 / ray.dir.x; // IEEE stability fix
	double divy = 1 / ray.dir.y;
	double divz = 1 / ray.dir.z;

	double tx0 = (octree->boxMin.x - ray.origin.x) * divx;
	double tx1 = (octree->boxMax.x - ray.origin.x) * divx;
	double ty0 = (octree->boxMin.y - ray.origin.y) * divy;
	double ty1 = (octree->boxMax.y - ray.origin.y) * divy;
	double tz0 = (octree->boxMin.z - ray.origin.z) * divz;
	double tz1 = (octree->boxMax.z - ray.origin.z) * divz;

	if (max(max(tx0, ty0), tz0) <= min(min(tx1, ty1), tz1))
		return true;
	else
		return false;
}


//   ray가 만나는 공간에 있는 오브젝트들을 골라낸다.
void RayTreeTraversal(TreeNode* node, Ray ray, vector<object*> *newObj)
{
	if (RayTraversal(node, ray))
	{
		int a = 0;

		for (int i = 0; i < 8; i++)
		{
			//   자식 node가 있다면 자식 node 먼저 visit
			if (node->children[i] != nullptr)
			{
				//cout << "hi" << endl;
				RayTreeTraversal(node->children[i], ray, newObj);
				a++;
			}
		}

		if (a == 0)
		{
			//cout << "hi" << endl;
			if (node->objs != nullptr)
			{
				//cout << newObj->size() << endl;
				//cout << "hello" << endl;
				vector<object*> ob = *(node->objs);

				//cout << "hi" << endl;
				for (int i = 0; i < ob.size(); i++)
				{
					newObj->push_back(ob[i]);
				}
			}
		}
	}

}

TreeNode *root = new TreeNode;


VECTOR3D raytrace(Ray ray, int depth)
{
	float min_t = 10000;
	object *o = nullptr;

	//Octree
	if (isoctree == 1)
	{
		vector<object*> newObjects;
		RayTreeTraversal(root, ray, &newObjects);

		for (int i = 0; i < newObjects.size(); i++)
		{
			float t;
			object *obj = newObjects[i];

			if (obj->hit(ray, &t))
			{
				if (t <= min_t)
				{
					//cout << "hello\n" << endl;
					min_t = t;
					o = obj;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < objects.size(); i++)
		{
			float t;
			object *obj = objects[i];

			if (obj->hit(ray, &t))
			{
				if (t <= min_t)
				{
					//cout << "hello\n" << endl;
					min_t = t;
					o = obj;
				}
			}
		}
	}

	if (min_t == 10000)
		return VECTOR3D(0.0f, 0.0f, 0.0f);

	//   shadow
	float shadow = 1.0;
	float min_st = 10000;
	VECTOR3D point = ray.origin + min_t * ray.dir;

	VECTOR3D L = light - point;
	L.Normalize();
	Ray shadow_ray(point, L);

	if (isoctree == 1)
	{
		vector<object*> newSobjects;
		RayTreeTraversal(root, shadow_ray, &newSobjects);

		for (int i = 0; i < newSobjects.size(); i++) {
			float t;
			object *obj = newSobjects[i];

			if (obj->hit(shadow_ray, &t)) {
				//cout << "dd" << endl;
				if (t <= min_st) {
					min_st = t;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < objects.size(); i++) {
			float t;
			object *obj = objects[i];

			if (obj->hit(shadow_ray, &t)) {
				//cout << "dd" << endl;
				if (t <= min_st) {
					min_st = t;
				}
			}
		}
	}

	if (min_st == 10000)
		shadow = 1.0;
	else
		shadow = 0.3;

	if (depth > 0)
	{
		VECTOR3D N = o->get_normal(point);
		N.Normalize();
		VECTOR3D Reflection = 2 * (N.InnerProduct(-1 * ray.dir))*N + ray.dir;
		Reflection.Normalize();

		int aa = (o->get_mode());

		if (aa != 1)
		{
			//cout << "s" << endl;
			return o->k_ambient + shadow * (o->getColor(point, light, ray.origin))
				+ 0.3*raytrace(Ray(point, Reflection), depth - 1);
			//+ 0.3*raytrace(Ray(point, Refraction), depth - 1);
		}
		else
		{
			return 1.0*raytrace(Ray(point, Reflection), depth - 1);
		}

	}
	else
		return o->k_ambient + shadow * o->getColor(point, light, ray.origin);
}


void ComputeRayDirection(int x, int y)
{
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	for (int i = 0; i < 800; i++)
	{
		for (int j = 0; j < 500; j++)
		{
			double nearX, nearY, nearZ;

			gluUnProject(i, j, 0, mvMatrix, projMatrix, viewport, &nearX, &nearY, &nearZ);

			VECTOR3D near((float)nearX, (float)nearY, (float)nearZ);

			Ray ray(eye, near - eye);
			ray.dir.Normalize();

			VECTOR3D color = raytrace(ray, depth);
			//cout << near.z << endl;

			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glLoadIdentity();

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();

			glColor3f(color.x, color.y, color.z);

			glBegin(GL_POINTS);
			//glVertex3d(near.x, near.y, near.z);
			glVertex3d(((float)i / 800.0f - 0.5f) * 2, ((float)j / 500.0f - 0.5f) * 2, 0.0f);
			glEnd();

			glMatrixMode(GL_MODELVIEW);
			glPopMatrix();
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		}
	}

	std::cout << "done" << endl;
}

void makeSphere(float d, float div) {

	for (int k = 0; k < div; k++) {
		for (int j = 0; j < div; j++) {
			for (int i = 0; i < div; i++) {
				
				//cout << -3.5*(d / div) + (d / div)*i << endl;

				sphere *tmp = new sphere();

				tmp->mode = 0;
				tmp->cen = { (float)((-d/2)+(d /(2*div)) + (d / div)*i),(float)((-d / 2) + (d / (2 * div)) + (d / div)*j),(float)((-d / 2) + (d / (2 * div)) + (d / div)*k) };
				tmp->rad = (d / 8)* 0.3;//개당 반지름 크기
				tmp->k_ambient = { 0.1f,0.1f,0.1f };
				tmp->k_diffuse = { 1.0f,1.0f,1.5f };
				tmp->k_specular = { 5.0f,5.0f,50.0f };
				tmp->k_shineness = 25.0f;

				objects.push_back(tmp);
			}
		}
	}

	cout << "생성" << endl;
}

void makerRoom(float rSize) {


	triangle *tmp = new triangle();

	tmp->mode = 1;
	tmp->a = { (rSize / 2),(rSize / 2),-(rSize / 2) };
	tmp->b = { -(rSize / 2),(rSize / 2),-(rSize / 2) };
	tmp->c = { -(rSize / 2),-(rSize / 2),-(rSize / 2) };

	tmp->va = { 0.0 , 0.0 , 1.0 };
	tmp->vb = { 0.0 , 0.0 , 1.0 };
	tmp->vc = { 0.0, 0.0 ,1.0 };

	tmp->k_ambient = { 0.1f,0.1f,0.2f };
	tmp->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp->k_specular = { 1.0f,1.0f,1.0f };
	tmp->k_shineness = 25.0f;

	objects.push_back(tmp);

	triangle *tmp1 = new triangle();

	tmp1->mode = 1;
	tmp1->a = { -(rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp1->b = { (rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp1->c = { (rSize / 2),(rSize / 2),-(rSize / 2) };

	tmp1->va = { 0.0 , 0.0 , 1.0 };
	tmp1->vb = { 0.0 , 0.0 , 1.0 };
	tmp1->vc = { 0.0, 0.0 ,1.0 };

	tmp1->k_ambient = { 0.1f,0.1f,0.2f };
	tmp1->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp1->k_specular = { 1.0f,1.0f,1.0f };
	tmp1->k_shineness = 25.0f;

	objects.push_back(tmp1);

	//우측면
	triangle *tmp2 = new triangle();

	tmp2->mode = 1;
	tmp2->a = { (rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp2->b = { (rSize / 2),(rSize / 2),(rSize / 2) };
	tmp2->c = { (rSize / 2),(rSize / 2),-(rSize / 2) };

	tmp2->va = { -1.0 , 0.0 , 0.0 };
	tmp2->vb = { -1.0 , 0.0 , 0.0 };
	tmp2->vc = { -1.0, 0.0 ,0.0 };

	tmp2->k_ambient = { 0.1f,0.1f,0.2f };
	tmp2->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp2->k_specular = { 1.0f,1.0f,1.0f };

	tmp2->k_shineness = 25.0f;

	objects.push_back(tmp2);

	triangle *tmp3 = new triangle();

	tmp3->mode = 1;
	tmp3->a = { (rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp3->b = { (rSize / 2),-(rSize / 2),(rSize / 2) };
	tmp3->c = { (rSize / 2),(rSize / 2),(rSize / 2) };

	tmp3->va = { -1.0 , 0.0 , 0.0 };
	tmp3->vb = { -1.0 , 0.0 , 0.0 };
	tmp3->vc = { -1.0, 0.0 ,0.0 };

	tmp3->k_ambient = { 0.1f,0.1f,0.2f };
	tmp3->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp3->k_specular = { 1.0f,1.0f,1.0f };
	tmp3->k_shineness = 25.0f;

	objects.push_back(tmp3);

	//왼쪽면

	triangle *tmp4 = new triangle();

	tmp4->mode = 1;
	tmp4->a = { -(rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp4->c = { -(rSize / 2),-(rSize / 2),(rSize / 2) };
	tmp4->b = { -(rSize / 2),(rSize / 2),(rSize / 2) };

	tmp4->va = { 1.0 , 0.0 , 0.0 };
	tmp4->vb = { 1.0 , 0.0 , 0.0 };
	tmp4->vc = { 1.0, 0.0 ,0.0 };

	tmp4->k_ambient = { 0.1f,0.1f,0.2f };
	tmp4->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp4->k_specular = { 1.0f,1.0f,1.0f };
	tmp4->k_shineness = 25.0f;

	objects.push_back(tmp4);

	triangle *tmp5 = new triangle();

	tmp5->mode = 1;
	tmp5->a = { -(rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp5->c = { -(rSize / 2),(rSize / 2),(rSize / 2) };
	tmp5->b = { -(rSize / 2),(rSize / 2),-(rSize / 2) };

	tmp5->va = { 1.0 , 0.0 , 0.0 };
	tmp5->vb = { 1.0 , 0.0 , 0.0 };
	tmp5->vc = { 1.0, 0.0 ,0.0 };

	tmp5->k_ambient = { 0.1f,0.1f,0.2f };
	tmp5->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp5->k_specular = { 1.0f,1.0f,1.0f };
	tmp5->k_shineness = 25.0f;

	objects.push_back(tmp5);

	//뒷면
	triangle *tmp6 = new triangle();

	tmp6->mode = 1;
	tmp6->a = { (rSize / 2),(rSize / 2),(rSize / 2) };
	tmp6->c = { -(rSize / 2),(rSize / 2),(rSize / 2) };
	tmp6->b = { -(rSize / 2),-(rSize / 2),(rSize / 2) };

	tmp6->va = { 0.0 , 0.0 , -1.0 };
	tmp6->vb = { 0.0 , 0.0 , -1.0 };
	tmp6->vc = { 0.0, 0.0 ,-1.0 };

	tmp6->k_ambient = { 0.1f,0.1f,0.2f };
	tmp6->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp6->k_specular = { 1.0f,1.0f,1.0f };
	tmp6->k_shineness = 25.0f;

	objects.push_back(tmp6);

	triangle *tmp7 = new triangle();

	tmp7->mode = 1;
	tmp7->a = { -(rSize / 2),-(rSize / 2),(rSize / 2) };
	tmp7->c = { (rSize / 2),-(rSize / 2),(rSize / 2) };
	tmp7->b = { (rSize / 2),(rSize / 2),(rSize / 2) };

	tmp7->va = { 0.0 , 0.0 , -1.0 };
	tmp7->vb = { 0.0 , 0.0 , -1.0 };
	tmp7->vc = { 0.0, 0.0 ,-1.0 };

	tmp7->k_ambient = { 0.1f,0.1f,0.2f };
	tmp7->k_diffuse = { 0.2f,0.2f,0.4f };
	tmp7->k_specular = { 1.0f,1.0f,1.0f };
	tmp7->k_shineness = 25.0f;

	objects.push_back(tmp7);


	//바닥
	triangle *tmp8 = new triangle();

	tmp8->mode = 0;
	tmp8->a = { (rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp8->b = { -(rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp8->c = { -(rSize / 2),-(rSize / 2),(rSize / 2) };

	tmp8->va = { 0.0 , 1.0 , 0.0 };
	tmp8->vb = { 0.0 , 1.0 , 0.0 };
	tmp8->vc = { 0.0, 1.0 ,0.0 };

	tmp8->k_ambient = { 0.1f,0.1f,0.1f };
	tmp8->k_diffuse = { 0.45f,0.45f,0.75f };
	tmp8->k_specular = { 1.0f,1.0f,1.5f };
	tmp8->k_shineness = 25.0f;

	objects.push_back(tmp8);

	triangle *tmp9 = new triangle();

	tmp9->mode = 0;
	tmp9->a = { (rSize / 2),-(rSize / 2),-(rSize / 2) };
	tmp9->b = { -(rSize / 2),-(rSize / 2),(rSize / 2) };
	tmp9->c = { (rSize / 2),-(rSize / 2),(rSize / 2) };

	tmp9->va = { 0.0 , 1.0 , 0.0 };
	tmp9->vb = { 0.0 , 1.0 , 0.0 };
	tmp9->vc = { 0.0, 1.0 ,0.0 };

	tmp9->k_ambient = { 0.1f,0.1f,0.1f };
	tmp9->k_diffuse = { 0.45f,0.45f,0.75f };
	tmp9->k_specular = { 1.0f,1.0f,1.5f };
	tmp9->k_shineness = 25.0f;

	objects.push_back(tmp9);


}


void Keyboard(unsigned char key, int x, int y)
{
	//To Do
	switch (key) {
		//카메라 이동
	case 'W':
	case 'w':
		if (cameraPos > 1)
			cameraPos -= 5.0f;
		cameraX = cameraPos * sin(cameraAngle*3.141592 / 180);
		cameraY = cameraPos / 3;
		cameraZ = cameraPos * cos(cameraAngle*3.141592 / 180);

		eye.x = cameraPos * sin(cameraAngle*3.141592 / 180);
		eye.y = cameraPos / 3;
		eye.z = cameraPos * cos(cameraAngle*3.141592 / 180);
		break;
	case 'S':
	case 's':
		cameraPos += 5.0f;
		cameraX = cameraPos * sin(cameraAngle*3.141592 / 180);
		cameraY = cameraPos / 3;
		cameraZ = cameraPos * cos(cameraAngle*3.141592 / 180);

		eye.x = cameraPos * sin(cameraAngle*3.141592 / 180);
		eye.y = cameraPos / 3;
		eye.z = cameraPos * cos(cameraAngle*3.141592 / 180);
		break;

		//카메라 회전
	case 'A':
	case 'a':
		cameraAngle -= 5.0f;
		cameraX = cameraPos * sin(cameraAngle*3.141592 / 180);
		cameraY = cameraPos / 3;
		cameraZ = cameraPos * cos(cameraAngle*3.141592 / 180);

		eye.x = cameraPos * sin(cameraAngle*3.141592 / 180);
		eye.y = cameraPos / 3;
		eye.z = cameraPos * cos(cameraAngle*3.141592 / 180);
		break;
	case 'D':
	case 'd':
		cameraAngle += 5.0f;
		cameraX = cameraPos * sin(cameraAngle*3.141592 / 180);
		cameraY = cameraPos / 3;
		cameraZ = cameraPos * cos(cameraAngle*3.141592 / 180);

		eye.x = cameraPos * sin(cameraAngle*3.141592 / 180);
		eye.y = cameraPos / 3;
		eye.z = cameraPos * cos(cameraAngle*3.141592 / 180);
		break;
	case 'g':
	case 'G':
		start = 1;
		break;
}

	glutPostRedisplay();
}


char* timeToString(struct tm *t) {
	static char s[20];

	sprintf(s, "%04d-%02d-%02d %02d:%02d:%02d",
		t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
		t->tm_hour, t->tm_min, t->tm_sec
	);
	return s;
}

void TimerFunc(int value)
{
	if (cameraAngle == 360)
	{
		finishTime = clock();
		cout << "\n-----Rotate Done-----" << endl;

		struct tm *t;
		time_t timer;

		timer = time(NULL);    // 현재 시각을 초 단위로 얻기
		t = localtime(&timer); // 초 단위의 시간을 분리하여 구조체에 넣기

		cout << "현재시간 : " << timeToString(t) << endl;

		cout << "걸린시간 : " << (double)(finishTime - startTime) / CLOCKS_PER_SEC << endl;

		start = 0;
		cameraAngle = 0;
	}
	//idle상태이면 delta를 0.001씩 증가
	if (start == 1 && cameraAngle < 360)
	{
		if (cameraAngle == 0)
		{
			cout << "-----Rotate Start-----" << endl;
			startTime = clock();

			struct tm *t;
			time_t timer;

			timer = time(NULL);    // 현재 시각을 초 단위로 얻기
			t = localtime(&timer); // 초 단위의 시간을 분리하여 구조체에 넣기

			cout << "현재시간 : " << timeToString(t) << endl;
		}

		cameraAngle += 40.0f;
		cameraX = cameraPos * sin(cameraAngle*3.141592 / 180);
		cameraZ = cameraPos * cos(cameraAngle*3.141592 / 180);

		eye.x = cameraPos * sin(cameraAngle*3.141592 / 180);
		eye.z = cameraPos * cos(cameraAngle*3.141592 / 180);


		glutPostRedisplay();
	}

	glutTimerFunc(100, TimerFunc, 1);
}

void Rendering(void)
{
	// 화면 버퍼 클리어
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// 화면을 제대로 바라보기 위해 카메라를 회전 후 이동
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX, cameraY, cameraZ, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	glPushMatrix();
	ComputeRayDirection(0, 0);
	glPopMatrix();

	// back 버퍼에 랜더링한 후 swap
	glutSwapBuffers();
}

void Reshape(int w, int h)
{
	// 뷰포트 설정
	glViewport(0, 0, w, h);

	// 원근 투영 사용
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, (float)w / h, 0.1, 500);

	// 모델뷰 매트릭스 초기화
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void EventHandlingAndLoop()
{
	glutKeyboardFunc(Keyboard);
	glutDisplayFunc(Rendering);  // 변환된 값에 따른 Rendering Callback 함수 등록
	glutReshapeFunc(Reshape);    // 윈도우 창 크기가 바뀌었을때 호출되는 Callback 함수 등록
	glutTimerFunc(1, TimerFunc, 1);

	glutMainLoop(); // 등록된 callback 함수를 반복하여 호출
}

int main(int argc, char** argv)
{
	//Meshload();
	//ComputeNormal();

	//RenderMesh(mesh3);

	isoctree = 1;   // 1이면 octree사용

					//Meshload();
					//ComputeNormal();
					//RenderMesh();
	makerRoom(200.0);


	makeSphere(40.0, 10.0);

	
	


	float size = 110.0f;
	root->boxMin.x = -1.0 * size;
	root->boxMin.y = -1.0 * size;
	root->boxMin.z = -1.0 * size;
	root->boxMax.x = size;
	root->boxMax.y = size;
	root->boxMax.z = size;

	if (isoctree == 1)
	{
		SpaceDivision(root, objects);
		//Subdivide(root, 1);
		cout << "Division complete" << endl;
	}

	//PostOrderTraversal(root);

	Initialize(argc, argv);      // 윈도우 생성, 배경색 설정

	EventHandlingAndLoop();      // Event Handling 및 Loop

								 // 에러 없이 끝났을 경우 0을 리턴함
	return 0;
}