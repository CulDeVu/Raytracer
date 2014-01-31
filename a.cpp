#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <limits>
#include <time.h>
#include <chrono>

using namespace std;

const float maxFloat = 80000000000000;
const int imageWidth = 255,
		  imageHeight = 255;

const int NUM_SAMPLES = 1,
		  NUM_BOUNCES = 1,
		  NUM_SCATTERS = 16;

enum MATERIAL_TYPE
{
	DIFFUSE, MIRROR
};

struct vec3
{
	float x, y, z;

	float length()
	{
		return sqrt(x*x + y*y + z*z);
	}

	vec3 normalized()
	{
		float len = length();
		vec3 me = { x/len, y/len, z/len };
		return me;
	}

	vec3 operator+(vec3 other) const
	{
		vec3 t = { other.x + x, other.y + y, other.z + z };
		return t;
	}
	vec3 operator-(vec3 other) const
	{
		vec3 t = { x - other.x, y - other.y, z - other.z };
		return t;
	}
	vec3 operator*(float f) const
	{
		vec3 t = { x * f, y * f, z * f };
		return t;
	}
};


vec3 newcolor(float r, float g, float b)
{
	vec3 temp = {r, g, b};
	if (r > 1 || g > 1 || b > 1 ||
		r < 0 || b < 0 || b < 0) 
	{
		cout << "SOMETHING WENT WRONG: (" << r << ", " << g << ", " << b << ")" << endl;
		exit(0);
	}
	return temp;
}

vec3 newray(float x, float y, float z)
{
	vec3 temp = {x, y, z};
	return temp;
}

float dot(vec3 v1, vec3 v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

float rdown(float f)
{
	if (f > 0)
		return (int)f;
	else
		return (int)f - 1;
}

class shape
{
public:
	virtual float intersect(vec3 o, vec3 ray) {}
	virtual vec3 getNorm(vec3 pt) {}
	virtual MATERIAL_TYPE getMaterialType() {}
	virtual vec3 getBRDF(vec3 pos) {}
	virtual vec3 getEmmision(vec3 pos) {}
};

class plane : public shape
{
public:
	plane(vec3 p, vec3 n)
	{
		pos = p;
		norm = n;
	}

	virtual float intersect(vec3 o, vec3 ray)
	{
		if (dot(ray, norm) == 0)
			return 0;
		return (dot(pos - o, norm) / dot(ray, norm));
	}

	virtual vec3 getNorm(vec3 pt)
	{
		return norm;
	}

	virtual MATERIAL_TYPE getMaterialType()
	{
		return DIFFUSE;
	}

	virtual vec3 getBRDF(vec3 pos)
	{
		vec3 color;
		if ((int)floor(pos.x) % 2 == 0 && (int)floor(pos.z) % 2 != 0)
			color = newcolor(1, 0, 0);
		else if ((int)floor(pos.x) % 2 != 0 && (int)floor(pos.z) % 2 == 0)
			color = newcolor(1, 0, 0);
		else
			color = newcolor(1, 1, 1);

		return color * (1 / 3.14159f);
	}

	virtual vec3 getEmmision(vec3 pos)
	{
		return newcolor(0, 0, 0);
	}

protected:
	vec3 pos, norm;
};

class circle : public shape
{
public:
	circle(vec3 p, float r)
	{
		pos = p;
		emm = newcolor(0, 0, 0);
		radius = r;
	}

	void setEmmision(vec3 color)
	{
		emm = color;
	}

	void setColor(vec3 c)
	{
		color = c;
	}

	void setMaterial(MATERIAL_TYPE t)
	{
		mat_type = t;
	}

	virtual float intersect(vec3 o, vec3 ray)
	{
		vec3 nray = ray.normalized();
		float a = dot(ray, ray);
		float b = 2 * dot(ray, o - pos);
		float c = dot(o - pos, o - pos) - radius*radius;
		if (b * b - 4 * a * c < 0)
			return 0;
		float t0 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a),
			  t1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
		return min(t0, t1);
	}

	virtual vec3 getNorm(vec3 pt)
	{
		vec3 norm = pt - pos;
		norm = norm.normalized();
		return norm;
	}

	virtual MATERIAL_TYPE getMaterialType()
	{
		return mat_type;
	}

	virtual vec3 getBRDF(vec3 pos)
	{
		return color * (1 / 3.14159f);
	}

	virtual vec3 getEmmision(vec3 pos)
	{
		return emm;
	}

protected:
	vec3 pos;
	vec3 emm;
	float radius;
	vec3 color;
	MATERIAL_TYPE mat_type;
};

vec3 buffer[imageWidth][imageHeight];
vector<shape*> shapes;

vec3 sample(vec3 o, vec3 ray, float bounces)
{
	float closestT = maxFloat;
	shape* closestShape = 0;

	for (int i = 0; i < shapes.size(); ++i)
	{
		float t = shapes[i]->intersect(o, ray);

		if (t > 0 && t < closestT)
		{
			closestT = t;
			closestShape = shapes[i];
		}
	}

	closestT -= 0.01f; // hackify :)

	if (closestShape == 0)
		return newcolor(0.0f, 0.0f, 0.0f);

	vec3 intersected = o + ray * closestT;
	vec3 norm = closestShape->getNorm(intersected);
	vec3 brdf = closestShape->getBRDF(intersected);

	vec3 emmisionTerm = closestShape->getEmmision(intersected);
	vec3 reflectTerm = newcolor(0, 0, 0);

	int numSamples = NUM_SCATTERS;
	if (closestShape->getMaterialType() == MIRROR)
		numSamples = 1;

	for (int i = 0; i < numSamples && bounces > 0; ++i)
	{
		vec3 castRay = newray(0, 0, 0);

		if (closestShape->getMaterialType() == MIRROR)
		{
			castRay = ray - norm * 2 * dot(ray, norm);
		}
		else
		{
			/*if (bounces == 1)
			{
				vec3 ptlight = newray(0, 5.0f, -11);
				vec3 d = ptlight - intersected;
				castRay = d;
			}
			else*/
			{
				while (true)
				{
					/*float rx = 2 * (float)rand() / RAND_MAX - 1.0f,
						  ry = 2 * (float)rand() / RAND_MAX - 1.0f,
						  rz = 2 * (float)rand() / RAND_MAX - 1.0f;
					if (rx*rx + ry*ry + rz*rz > 1.0f)
						continue;
					castRay = newray(rx, ry, rz);*/

					// cylinder-to-sphere projection of random stuff
					float h = 2 * (float)rand() / RAND_MAX - 1.0f;
					float angle = 2 * 3.14159 * (float)rand() / RAND_MAX;
					float r = sqrt(1 - h*h);

					float rz = h,
						  rx = r * cos(angle),
						  ry = r * sin(angle);

					castRay = newray(rx, ry, rz);

					if (dot(castRay, norm) >= 0)
						break;
				}
			}
		}

		float cosAngle = dot(castRay.normalized(), norm);
		if (cosAngle < 0)
			cosAngle = 0;

		vec3 temp = sample(intersected, castRay, bounces - 1);

		reflectTerm = reflectTerm + newray(brdf.x*temp.x, brdf.y*temp.y, brdf.z*temp.z) * cosAngle;
	}

	reflectTerm = reflectTerm * (1.0f /  (2 * 3.14159f * numSamples));

	vec3 color = emmisionTerm + reflectTerm;

	return color;
}

int main() {
	for (int y = 0; y < imageHeight; ++y)
		for (int x = 0; x < imageWidth; ++x)
			buffer[x][y] = newcolor(0.0f, 0.0f, 0.0f);

	plane* bottom = new plane(newray(0, -3, 0), newray(0, 1, 0));
	shapes.push_back(bottom);

	plane* p1 = new plane(newray(0, 0, -20), newray(0, 0, 1));
	shapes.push_back(p1);

	plane* p2 = new plane(newray(-5, 0, 0), newray(1, 0, 0));
	shapes.push_back(p2);

	plane* right = new plane(newray(5, 0, 0), newray(-1, 0, 0));
	shapes.push_back(right);

	plane* top = new plane(newray(0, 5, 0), newray(0, -1, 0));
	shapes.push_back(top);

	plane* back = new plane(newray(0, 0, 5), newray(0, 0, -1));
	shapes.push_back(back);

	circle* ball1 = new circle(newray(-3, 0, -15), 1.5f);
	ball1->setColor(newcolor(1, 1, 1));
	ball1->setMaterial(DIFFUSE);
	shapes.push_back(ball1);

	circle* light = new circle(newray(0, 12.8f, -11), 8.0f);
	light->setColor(newcolor(0, 0, 0));
	light->setEmmision(newray(1000, 1000, 1000));
	shapes.push_back(light);

	chrono::time_point<chrono::system_clock> start, end;
	start = chrono::system_clock::now();

	//#pragma omp parallel for schedule(dynamic)
	for (int y = -imageHeight/2; y < imageHeight/2; ++y)
	{
		//if (y == 0)
		//	cout << "50% done" << endl;
		for (float x = -imageWidth/2; x < imageWidth/2; ++x)
		{
			vec3 totalColor = newcolor(0, 0, 0);

			for (int i = 0; i < NUM_SAMPLES; ++i)
			{
				/*float rx = 0.5f * ((float)rand() / RAND_MAX - 0.5f) / (imageWidth),
					  ry = 0.5f * ((float)rand() / RAND_MAX - 0.5f) / (imageHeight);*/
				float rx = 0, ry = 0;
				vec3 ray = newray(x/2.0f/imageWidth + rx, y/2.0f/imageHeight + ry, -0.5f);

				vec3 color = sample(newray(0, 0, 0), ray, NUM_BOUNCES);

				totalColor = totalColor + color;
			}

			if (totalColor.x < 0 || totalColor.y < 0 || totalColor.z < 0)
				cout << "adfs;jkl" << endl;
			totalColor = totalColor * (1.0f / NUM_SAMPLES);
			if (totalColor.x > 1.0f)
				totalColor.x = 1.0f; //totalColor.x / (totalColor.x + 1);
			if (totalColor.y > 1.0f)
				totalColor.y = 1.0f; //totalColor.y / (totalColor.y + 1);
			if (totalColor.z > 1.0f)
				totalColor.z = 1.0f; //totalColor.z / (totalColor.z + 1);
			buffer[(int)(x + imageWidth/2)][(int)(y + imageHeight/2)] = totalColor;
		}
	}

	end = chrono::system_clock::now();

	chrono::duration<double> elapsed_seconds = end - start;
	cout << "Elapsed time: " << elapsed_seconds.count() << endl;
	cout << "Time per pixel: " << 1000 * elapsed_seconds.count() / (imageWidth*imageHeight) << endl;
	//printf("Elapsed time: %.2lf seconds\n", elapsed_seconds.count());

	ofstream file("image2.ppm");
	file << "P3 " << imageWidth << " " << imageHeight << " 255" << endl;


	for (int y = imageHeight - 1; y >= 0; --y)
		for (int x = 0; x < imageWidth; ++x)
			file << (int)(buffer[x][y].x*255) << " " << (int)(buffer[x][y].y*255) << " " << (int)(buffer[x][y].z*255) << " ";

	file.close();

	return 0;
}
