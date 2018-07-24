#include "Mesh.h"

using namespace std;

void Mesh::LoadMesh(string str)
{
	ifstream file(str);
	string inputString;
	while (!file.eof())
	{
		
		getline(file, inputString);
		if (inputString[0] == '#')
			continue;
		string delimiter = " ";
		string tokens[4];
		size_t pos = 0;
		int index = 0;
		string token;
		while ((pos = inputString.find(delimiter)) != std::string::npos)
		{
			
			token = inputString.substr(0, pos);
			tokens[index] = token;
			inputString.erase(0, pos + delimiter.length());
			index++;
		}
		tokens[3] = inputString;
		if (tokens[0] == "v")
		{
			Vertex v;
			v.position = VECTOR3D(stof(tokens[1]), stof(tokens[2]), stof(tokens[3]));
			vertexArray.push_back(v);
		}
		else if (tokens[0] == "f")
		{
			Face f;
			f.vertex0 = stoi(tokens[1]);
			f.vertex1 = stoi(tokens[2]);
			f.vertex2 = stoi(tokens[3]);
			f.vertex0--;
			f.vertex1--;
			f.vertex2--;
			faceArray.push_back(f);
		}
	}
	file.close();
	
}

void Mesh::ComputeFaceNormal()
{
	for (int i = 0; i < faceArray.size(); i++)
	{
		Face &f = faceArray[i];

		VECTOR3D v0 = vertexArray[f.vertex0].position;
		VECTOR3D v1 = vertexArray[f.vertex1].position;
		VECTOR3D v2 = vertexArray[f.vertex2].position;

		VECTOR3D va = v1 - v0;
		VECTOR3D vb = v2 - v0;
		VECTOR3D vc = va.CrossProduct(vb);
		VECTOR3D vd = va.CrossProduct(vb);
		vc.Normalize();
		f.normal = vc;
	}
}

void Mesh::FindNeighborFaceArray()
{
	for (int i = 0; i < vertexArray.size(); i++)
	{
		for (int j = 0; j < faceArray.size(); j++)
		{
			if (faceArray[j].vertex0 == i || faceArray[j].vertex1 == i || faceArray[j].vertex2 == i)
			{
				vertexArray[i].neighborFaces.push_back(j);
			}
		}
	}
}

void Mesh::ComputeVertexNormal()
{
	for (int i = 0; i < vertexArray.size(); i++)
	{
		Vertex& v = vertexArray[i];
		for (int j = 0; j < vertexArray[i].neighborFaces.size(); j++)
		{
			v.normal += faceArray[vertexArray[i].neighborFaces[j]].normal;
		}
		v.normal /= vertexArray[i].neighborFaces.size();
	}
}