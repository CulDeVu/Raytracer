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
const int imageWidth = 728,
		  imageHeight = 728;

const int NUM_SAMPLES = 16,
		  NUM_BOUNCES = 3,
		  NUM_SCATTERS = 1;

enum Matterial
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

class shape
{
public:
	virtual float intersect(vec3 o, vec3 ray) {}
	virtual vec3 getNorm(vec3 pt) {}
	virtual vec3 getBRDF(vec3 pos) {}
	virtual Matterial getMat() {return DIFFUSE;}
	virtual vec3 getEmmision(vec3 pos) {}
};

class plane : public shape
{
public:
	plane(vec3 p, vec3 n, vec3 c)
	{
		pos = p;
		norm = n;
		color = c;
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

	virtual vec3 getBRDF(vec3 pos)
	{
		/*vec3 color;
		if ((int)floor(pos.x) % 2 == 0 && (int)floor(pos.z) % 2 != 0)
			color = newcolor(1, 0, 0);
		else if ((int)floor(pos.x) % 2 != 0 && (int)floor(pos.z) % 2 == 0)
			color = newcolor(1, 0, 0);
		else
			color = newcolor(1, 1, 1);*/

		return color * (1 / 3.14159f);
	}

	virtual vec3 getEmmision(vec3 pos)
	{
		return newcolor(0, 0, 0);
	}

protected:
	vec3 pos, norm;
	vec3 color;
};

class circle : public shape
{
public:
	circle(vec3 p, float r)
	{
		pos = p;
		emm = newcolor(0, 0, 0);
		radius = r;
		mat = DIFFUSE;
	}

	void setEmmision(vec3 color)
	{
		emm = color;
	}

	void setColor(vec3 c)
	{
		color = c;
	}

	void setMat(Matterial m)
	{
		mat = m;
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

	virtual Matterial getMat()
	{
		return mat;
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
	Matterial mat;
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

	// i hate floating point numbers >:(
	int floatbits = *(int*)&closestT - 1;
	closestT = *(float*)&(floatbits);

	if (closestShape == 0)
		return newcolor(0.0f, 0.0f, 0.0f);

	vec3 intersected = o + ray * closestT;
	vec3 norm = closestShape->getNorm(intersected);
	vec3 brdf = closestShape->getBRDF(intersected);

	vec3 emmisionTerm = closestShape->getEmmision(intersected);
	vec3 reflectTerm = newcolor(0, 0, 0);

	float numSamples = 1;
	float PDF = 1;

	if (closestShape->getMat() == MIRROR)
	{
		vec3 castRay = ray - norm * 2 * dot(ray, norm);

		float cosAngle = dot(castRay.normalized(), norm);
		if (cosAngle < 0)
			cosAngle = 0;

		vec3 temp = sample(intersected, castRay, bounces - 1);

		reflectTerm = reflectTerm + newray(brdf.x*temp.x, brdf.y*temp.y, brdf.z*temp.z) * cosAngle;
	}
	else
	{
		numSamples = NUM_SCATTERS;
		PDF = 1 / (2 * 3.14159f);
		for (int i = 0; i < numSamples && bounces > 0; ++i)
		{
			vec3 castRay = newray(0, 0, 0);

			while (true)
			{
				// technically, this is faster, sooooo :P
				float rx = 2 * (float)rand() / RAND_MAX - 1.0f,
					  ry = 2 * (float)rand() / RAND_MAX - 1.0f,
					  rz = 2 * (float)rand() / RAND_MAX - 1.0f;
				if (rx*rx + ry*ry + rz*rz > 1.0f)
					continue;
				castRay = newray(rx, ry, rz);

				if (dot(castRay, norm) >= 0)
					break;
			}

			float cosAngle = dot(castRay.normalized(), norm);
			if (cosAngle < 0)
				cosAngle = 0;

			vec3 temp = sample(intersected, castRay, bounces - 1);

			reflectTerm = reflectTerm + newray(brdf.x*temp.x, brdf.y*temp.y, brdf.z*temp.z) * cosAngle;
		}
	}

	reflectTerm = reflectTerm * (1.0f /  (PDF * numSamples));

	vec3 color = emmisionTerm + reflectTerm;

	return color;
}

int main() {
	for (int y = 0; y < imageHeight; ++y)
		for (int x = 0; x < imageWidth; ++x)
			buffer[x][y] = newcolor(0.0f, 0.0f, 0.0f);

	plane* bottom = new plane(newray(0, -3, 0), newray(0, 1, 0), newcolor(1, 1, 1));
	shapes.push_back(bottom);

	plane* front = new plane(newray(0, 0, -20), newray(0, 0, 1), newcolor(1, 1, 1));
	shapes.push_back(front);

	plane* left = new plane(newray(-5, 0, 0), newray(1, 0, 0), newcolor(1, 0.2f, 0.2f));
	shapes.push_back(left);

	plane* right = new plane(newray(5, 0, 0), newray(-1, 0, 0), newcolor(0.2f, 1, 0.2f));
	shapes.push_back(right);

	plane* top = new plane(newray(0, 5, 0), newray(0, -1, 0), newcolor(1, 1, 1));
	shapes.push_back(top);

	plane* back = new plane(newray(0, 0, 5), newray(0, 0, -1), newcolor(0.5f, 0.5f, 0.5f));
	shapes.push_back(back);

	circle* ball1 = new circle(newray(-3, 0, -15), 1.5f);
	ball1->setColor(newcolor(1, 1, 1));
	shapes.push_back(ball1);

	circle* ball2 = new circle(newray(3, -1, -10), 1.5f);
	ball2->setColor(newcolor(1, 1, 1));
	shapes.push_back(ball2);

	circle* ball3 = new circle(newray(0, -3, -20), 3.0f);
	ball3->setColor(newcolor(1, 1, 1));
	ball3->setMat(MIRROR);
	shapes.push_back(ball3);

	circle* light = new circle(newray(0, 12.8f, -11), 8.0f);
	light->setColor(newcolor(0, 0, 0));
	light->setEmmision(newray(15, 15, 15));
	shapes.push_back(light);

	auto start = chrono::high_resolution_clock::now();

	//#pragma omp parallel for schedule(dynamic)
	for (int y = -imageHeight/2; y < imageHeight/2; ++y)
	{
		for (float x = -imageWidth/2; x < imageWidth/2; ++x)
		{
			vec3 totalColor = newcolor(0, 0, 0);

			for (int i = 0; i < NUM_SAMPLES; ++i)
			{
				float rx = 1.0f * ((float)rand() / RAND_MAX - 0.5f) / (imageWidth),
					  ry = 1.0f * ((float)rand() / RAND_MAX - 0.5f) / (imageHeight);
				//float rx = 0, ry = 0;
				vec3 ray = newray(x/2.0f/imageWidth + rx, y/2.0f/imageHeight + ry, -0.5f);

				vec3 color = sample(newray(0, 0, 0), ray, NUM_BOUNCES);

				totalColor = totalColor + color;
			}

			if (totalColor.x < 0 || totalColor.y < 0 || totalColor.z < 0)
				cout << "adfs;jkl" << endl;
			totalColor = totalColor * (1.0f / NUM_SAMPLES);

			totalColor.x = totalColor.x / (totalColor.x + 1);
			totalColor.y = totalColor.y / (totalColor.y + 1);
			totalColor.z = totalColor.z / (totalColor.z + 1);

			buffer[(int)(x + imageWidth/2)][(int)(y + imageHeight/2)] = totalColor;
		}
	}

	auto end = chrono::high_resolution_clock::now();

	cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(end - start).count() << endl;
	cout << "Milliseconds per pixel: " <<  (float)chrono::duration_cast<chrono::milliseconds>(end - start).count() / (imageWidth*imageHeight) << endl;
	//printf("Elapsed time: %.2lf seconds\n", elapsed_seconds.count());

	ofstream file("image2.ppm");
	file << "P3 " << imageWidth << " " << imageHeight << " 255" << endl;


	for (int y = imageHeight - 1; y >= 0; --y)
		for (int x = 0; x < imageWidth; ++x)
			file << (int)(buffer[x][y].x*255) << " " << (int)(buffer[x][y].y*255) << " " << (int)(buffer[x][y].z*255) << " ";

	file.close();

	return 0;
}
