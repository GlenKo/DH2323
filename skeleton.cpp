// Done by:
// Glen Ko - glen_ko@hotmail.com
// Ian Leow - ianleow.nus@gmail.com

#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

#define PI 3.1415926535f

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::ivec2;
using glm::mat3;

// ----------------------------------------------------------------------------
// STRUCT DECLARATIONS

struct Pixel {
	int x;
	int y;
	float zinv;
};

struct PixelI {
	int x;
	int y;
	float zinv;
	vec3 illumination;
};

struct PixelP {
	int x;
	int y;
	float zinv;
	vec3 pos3d;
};

struct Vertex {
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};

Pixel operator-(Pixel const& lhs, Pixel const& rhs);
PixelI operator-(PixelI const& lhs, PixelI const& rhs);
PixelP operator-(PixelP const& lhs, PixelP const& rhs);

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;

float focalLength = 500.0f;
vec3 cameraPosition(-0.0065f, -0.1473f, -0.2671f);
float yaw = 0.0436f;
float pitch = 0.1614f;
float roll = 0.0f; //TODO

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

vec3 lightPos(0.1565f, -0.18f, 1.7829f);
//vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 14.f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.05f*vec3(1, 1, 1);

vec3 currentNormal;
vec3 currentReflectance;

// Loader variables
const int PLY_LOADER_OK = 0;
const int PLY_LOADER_INVALID_IDENTIFIER = 1;
const int PLY_LOADER_UNKNOWN_FORMAT = 2;
const int PLY_LOADER_INVALID_FORMATTING = 3;
const int PLY_LOADER_UNSUPPORTED_FACE_PROPERTY = 4;
const int PLY_LOADER_UNSUPPORTED_POLYGON = 5;
const int PLY_LOADER_WRONG_NUMBER_OF_AXIS = 6;

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader(const vec3& v, ivec2& p);

void DrawEdges();
void DrawPolygonEdges(const vector<vec3>& vertices);
void DrawLineSDL(ivec2 a, ivec2 b, vec3 color);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);

mat3 getRotationMatrix(bool negateYaw = false, bool negatePitch = false);

void DrawColoredPolygons();
void DrawPolygonWithColor(const vector<vec3>& vertices, vec3& color);
void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels);
void DrawColoredRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3& color);

void DrawColoredPolygonsWithDepth();
void DrawPolygonWithColorAndDepth(const vector<vec3>& vertices, vec3& color);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawColoredRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3& color);
void VertexShader(const vec3& v, Pixel& p);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);

void DrawWithIllumination();
void DrawPolygonWithIllumination(const vector<Vertex>& vertices);
void ComputePolygonRows(const vector<PixelI>& vertexPixels, vector<PixelI>& leftPixels, vector<PixelI>& rightPixels);
void DrawColoredRows(const vector<PixelI>& leftPixels, const vector<PixelI>& rightPixels);
void VertexShader(const Vertex& v, PixelI& p);
void PixelShader(const PixelI& p);
void Interpolate(PixelI a, PixelI b, vector<PixelI>& result);

void DrawWithPixelIllumination();
void DrawPolygonWithPixelIllumination(const vector<Vertex>& vertices);
void ComputePolygonRows(const vector<PixelP>& vertexPixels, vector<PixelP>& leftPixels, vector<PixelP>& rightPixels);
void DrawColoredRows(const vector<PixelP>& leftPixels, const vector<PixelP>& rightPixels);
void VertexShader(const Vertex& v, PixelP& p);
void PixelShader(const PixelP& p);
void Interpolate(PixelP a, PixelP b, vector<PixelP>& result);

template<class T> T myMax(T a, T b);
template<class T> T myMin(T a, T b);
PixelP myAbs(PixelP const& p);
PixelI myAbs(PixelI const& p);
Pixel myAbs(Pixel const& p);
float myAbs(float const& v);
int myAbs(int const& v);

// Functions Coded by Glen
bool isBackface(vec3 pos3d, vec3 triangleNormal);
void writeNormal(Triangle& triangle);

// ----------------------------------------------------------------------------
void LoadCustomModel(vector<Triangle>& triangles);

int LoadPlyFile(const char * filePath, vector<Triangle>& triangles, bool reverseX = false, bool reverseY = false, bool reverseZ = false);
char * readHeaderLine(FILE * input);
int systemIsBigEndian();
void reverseEndianness(unsigned char * binaryData, int size);
void printBounds(vector<Triangle>& triangles);

int main(int argc, char* argv[]) {
	//LoadCustomModel(triangles);
	
	int plyLoaderErrorCode = LoadPlyFile("happy_vrip.ply", triangles, false, true, true);
	if (plyLoaderErrorCode != PLY_LOADER_OK) {
		cout << "Unable to load file, error code " << plyLoaderErrorCode << endl;
		return 1;
	}

	printBounds(triangles);

	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	// ======================================== TEST 1 ========================================
	//while( NoQuitMessageSDL() ) {
	//	Update();
	//	Draw();
	//}
	//SDL_SaveBMP( screen, "test1.bmp" );

	// ======================================== TEST 2 ========================================
	//while (NoQuitMessageSDL()) {
	//	Update();
	//	DrawEdges();
	//}
	//SDL_SaveBMP(screen, "test2.bmp");

	// ======================================== TEST 3 ========================================
	//vector<ivec2> vertexPixels(3);
	//vertexPixels[0] = ivec2(10, 5);
	//vertexPixels[1] = ivec2(5, 10);
	//vertexPixels[2] = ivec2(15, 15);
	//vector<ivec2> leftPixels;
	//vector<ivec2> rightPixels;
	//ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	//for (int row = 0; row<leftPixels.size(); ++row)
	//{
	//	cout << "Start: ("
	//		<< leftPixels[row].x << ","
	//		<< leftPixels[row].y << "). "
	//		<< "End: ("
	//		<< rightPixels[row].x << ","
	//		<< rightPixels[row].y << "). " << endl;
	//}
	//system("pause");

	// ======================================== TEST 4 ========================================
	//while (NoQuitMessageSDL()) {
	//	Update();
	//	DrawColoredPolygons();
	//}
	//SDL_SaveBMP(screen, "test4.bmp");

	// ======================================== TEST 5 ========================================
	//while (NoQuitMessageSDL()) {
	//	Update();
	//	DrawColoredPolygonsWithDepth();
	//}
	//SDL_SaveBMP(screen, "test5.bmp");

	// ======================================== TEST 6 ========================================
	//while (NoQuitMessageSDL()) {
	//	Update();
	//	DrawWithIllumination();
	//}
	//SDL_SaveBMP(screen, "test6.bmp");

	// ======================================== TEST 7/8 ========================================
	while (NoQuitMessageSDL()) {
		Update();
		DrawWithPixelIllumination();
	}
	SDL_SaveBMP(screen, "test7.bmp");

	return 0;
}

void Update() {
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;
	//cout << lightPos.x << " " << lightPos.y << " " << lightPos.z << endl;

	Uint8* keystate = SDL_GetKeyState(0);


	bool shiftHeld = false;
	if (keystate[SDLK_RSHIFT] || keystate[SDLK_LSHIFT]) {
		shiftHeld = true;
	}

	bool controlHeld = false;
	if (keystate[SDLK_RCTRL] || keystate[SDLK_LCTRL]) {
		controlHeld = true;
	}

	if (keystate[SDLK_UP]) {
		if (shiftHeld && !controlHeld) {
			cameraPosition = cameraPosition - (vec3(0.0f, 0.01f, 0.0f) * getRotationMatrix(true, true));
		}
		if (!shiftHeld && controlHeld) {
			cameraPosition = cameraPosition + (vec3(0.0f, 0.0f, 0.01f) * getRotationMatrix(true, true));
		}
		if (!shiftHeld && !controlHeld) {
			pitch -= PI / 720.0f; // 0.25 degrees
		}
	}

	if (keystate[SDLK_DOWN]) {
		if (shiftHeld && !controlHeld) {
			cameraPosition = cameraPosition + (vec3(0.0f, 0.01f, 0.0f) * getRotationMatrix(true, true));
		}
		if (!shiftHeld && controlHeld) {
			cameraPosition = cameraPosition - (vec3(0.0f, 0.0f, 0.01f) * getRotationMatrix(true, true));
		}
		if (!shiftHeld && !controlHeld) {
			pitch += PI / 720.0f; // 0.25 degrees
		}
	}

	if (keystate[SDLK_RIGHT]) {
		if (shiftHeld) {
			cameraPosition = cameraPosition + (vec3(0.01f, 0.0f, 0.0f) * getRotationMatrix(true, true));
		}
		else {
			yaw -= PI / 720.0f; // 0.25 degrees
		}
	}

	if (keystate[SDLK_LEFT]) {
		if (shiftHeld) {
			cameraPosition = cameraPosition - (vec3(0.01f, 0.0f, 0.0f) * getRotationMatrix(true, true));
		}
		else {
			yaw += PI / 720.0f; // 0.25 degrees
		}
	}

	if (keystate[SDLK_w]) {
		lightPos = lightPos + vec3(0.0f, 0.0f, 0.05f);
	}

	if (keystate[SDLK_s]) {
		lightPos = lightPos - vec3(0.0f, 0.0f, 0.05f);
	}

	if (keystate[SDLK_d]) {
		lightPos = lightPos + vec3(0.05f, 0.0f, 0.0f);
	}

	if (keystate[SDLK_a]) {
		lightPos = lightPos - vec3(0.05f, 0.0f, 0.0f);
	}

	if (keystate[SDLK_e]) {}

	if (keystate[SDLK_q]) {}
}

void Draw() {
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (size_t i = 0; i<triangles.size(); ++i) {
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		for (int v = 0; v<3; ++v) {
			ivec2 projPos;
			VertexShader(vertices[v], projPos);
			vec3 color(1, 1, 1);
			PutPixelSDL(screen, projPos.x, projPos.y, color);
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void VertexShader(const vec3& v, ivec2& p) {
	vec3 vprime = v - cameraPosition;
	vprime = vprime * getRotationMatrix();
	p.x = int(focalLength * vprime.x / vprime.z + SCREEN_WIDTH / 2.0f);
	p.y = int(focalLength * vprime.y / vprime.z + SCREEN_HEIGHT / 2.0f);
}

void DrawEdges() {
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (size_t i = 0; i<triangles.size(); ++i) {
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		DrawPolygonEdges(vertices);
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void DrawPolygonEdges(const vector<vec3>& vertices) {
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<ivec2> projectedVertices(V);
	for (int i = 0; i<V; ++i) {
		VertexShader(vertices[i], projectedVertices[i]);
	}
	// Loop over all vertices and draw the edge from it to the next vertex:
	for (int i = 0; i<V; ++i) {
		int j = (i + 1) % V; // The next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(projectedVertices[i], projectedVertices[j], color);
	}
}

void DrawLineSDL(ivec2 a, ivec2 b, vec3 color) {
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<ivec2> result(pixels);
	Interpolate(a, b, result);

	for (size_t i = 0; i<result.size(); ++i) {
		PutPixelSDL(screen, result[i].x, result[i].y, color);
	}
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) {
	//int tt1 = SDL_GetTicks();

	int N = result.size();
	vec2 step = vec2(b - a) / myMax(N - 1.0f, 1.0f);
	vec2 current(a);
	for (int i = 0; i<N; ++i) {
		result[i] = current;
		current += step;
	}

	//TODO set up better
	// Set the last result to b, since it should be, but might not be due to float error calculations
	result[N - 1] = b;

	//cout << "debug" << endl
	//	<< "0: " << result[0].y << endl
	//	<< "N-1: " << result[N - 1].y << endl
	//	<< "first: " << a.y << endl
	//	<< "last: " << b.y << endl
	//	<< "step: " << step.y << endl
	//	<< "N: " << N << endl
	//	<< endl; cout.flush();

	//int tt2 = SDL_GetTicks();
	//int diff = tt2 - tt1;
	//if (diff >= 10) {
	//	cout << N << " ---------- " << a.x << " " << a.y << " ========== " << b.x << " " << b.y << " ++++++++++ " << step.x << " " << step.y << endl;
	//}
}

mat3 getRotationMatrix(bool negateYaw, bool negatePitch) {
	// row-major i.e. vec3 * rotationMatrix, as per computer graphics 'convention'

	float thisYaw = (negateYaw) ? -yaw : yaw;
	float thisPitch = (negatePitch) ? -pitch : pitch;

	return mat3(
		vec3(cos(thisYaw), 0.0f, sin(thisYaw)),
		vec3(0.0f, 1.0f, 0.0f),
		vec3(-sin(thisYaw), 0.0f, cos(thisYaw))
		) *
		mat3(
		vec3(1.0f, 0.0f, 0.0f),
		vec3(0.0f, cos(thisPitch), -sin(thisPitch)),
		vec3(0.0f, sin(thisPitch), cos(thisPitch))
		);
}

void DrawColoredPolygons() {
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (size_t i = 0; i<triangles.size(); ++i) {
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		DrawPolygonWithColor(vertices, triangles[i].color);
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

//aka DrawPolygon(...) mentioned in the lab instructions
void DrawPolygonWithColor(const vector<vec3>& vertices, vec3& color) {
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<ivec2> projectedVertices(V);
	for (int i = 0; i<V; ++i) {
		VertexShader(vertices[i], projectedVertices[i]);
	}

	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
	DrawColoredRows(leftPixels, rightPixels, color);
}

void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels) {
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();
	for (size_t i = 0; i < vertexPixels.size(); i++) {
		minY = myMin(minY, vertexPixels[i].y);
		maxY = myMax(maxY, vertexPixels[i].y);
	}
	int numOfRows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels = vector<ivec2>(numOfRows);
	rightPixels = vector<ivec2>(numOfRows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < numOfRows; i++) {
		leftPixels[i].x = numeric_limits<int>::max();
		leftPixels[i].y = minY + i;

		rightPixels[i].x = numeric_limits<int>::min();
		rightPixels[i].y = minY + i;
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (size_t i = 0; i < vertexPixels.size(); ++i) {
		int j = (i + 1) % vertexPixels.size(); // The next vertex

		ivec2 delta = glm::abs(vertexPixels[i] - vertexPixels[j]);
		int pixels = glm::max(delta.x, delta.y) + 1;

		vector<ivec2> result(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], result);
		for (size_t k = 0; k < result.size(); k++) {
			leftPixels[result[k].y - minY].x = myMin(leftPixels[result[k].y - minY].x, result[k].x);
			rightPixels[result[k].y - minY].x = myMax(rightPixels[result[k].y - minY].x, result[k].x);
		}
	}
}

void DrawColoredRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3& color) {
	for (size_t i = 0; i < leftPixels.size(); i++) {
		for (int j = leftPixels[i].x; j <= rightPixels[i].x; j++) {
			PutPixelSDL(screen, j, leftPixels[i].y, color);
		}
	}
}

// For test 5 ==============================================================================================

void DrawColoredPolygonsWithDepth() {
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	// Reset depthBuffer
	for (int i = 0; i < SCREEN_HEIGHT; i++)
		for (int j = 0; j < SCREEN_WIDTH; j++)
			depthBuffer[i][j] = 0.0f;

	for (size_t i = 0; i<triangles.size(); ++i) {
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		DrawPolygonWithColorAndDepth(vertices, triangles[i].color);
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void DrawPolygonWithColorAndDepth(const vector<vec3>& vertices, vec3& color) {
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<Pixel> projectedVertices(V);
	for (int i = 0; i < V; ++i) {
		VertexShader(vertices[i], projectedVertices[i]);
	}

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
	DrawColoredRows(leftPixels, rightPixels, color);
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();
	for (size_t i = 0; i < vertexPixels.size(); i++) {
		minY = myMin(minY, vertexPixels[i].y);
		maxY = myMax(maxY, vertexPixels[i].y);
	}
	int numOfRows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels = vector<Pixel>(numOfRows);
	rightPixels = vector<Pixel>(numOfRows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < numOfRows; i++) {
		leftPixels[i].x = numeric_limits<int>::max();
		leftPixels[i].y = minY + i;
		leftPixels[i].zinv = 0.0f;

		rightPixels[i].x = numeric_limits<int>::min();
		rightPixels[i].y = minY + i;
		leftPixels[i].zinv = 0.0f;
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (size_t i = 0; i < vertexPixels.size(); ++i) {
		int j = (i + 1) % vertexPixels.size(); // The next vertex

		Pixel delta = myAbs(vertexPixels[i] - vertexPixels[j]);
		int pixels = myMax(delta.x, delta.y) + 1;

		vector<Pixel> result(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], result);
		for (size_t k = 0; k < result.size(); k++) {
			int indexY = result[k].y - minY;
			if (leftPixels[indexY].x > result[k].x) {
				leftPixels[indexY] = result[k];
			}
			if (rightPixels[indexY].x < result[k].x) {
				rightPixels[indexY] = result[k];
			}
		}
	}
}

void DrawColoredRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3& color) {
	for (size_t i = 0; i < leftPixels.size(); i++) {
		int rowLength = rightPixels[i].x - leftPixels[i].x + 1;

		vector<Pixel> result(rowLength);
		Interpolate(leftPixels[i], rightPixels[i], result);

		for (int j = 0; j < rowLength; j++) {
			int x = result[j].x;
			int y = result[j].y;

			if (x < 0 || y < 0 || x >= SCREEN_WIDTH || y >= SCREEN_HEIGHT)
				continue;

			if (result[j].zinv >= depthBuffer[x][y]) {
				depthBuffer[x][y] = result[j].zinv;
				PutPixelSDL(screen, x, y, color);
			}
		}
	}
}

void VertexShader(const vec3& v, Pixel& p) {
	vec3 vprime = v - cameraPosition;
	vprime = vprime * getRotationMatrix();

	p.x = int(focalLength * vprime.x / vprime.z + SCREEN_WIDTH / 2.0f);
	p.y = int(focalLength * vprime.y / vprime.z + SCREEN_HEIGHT / 2.0f);
	p.zinv = 1.0f / vprime.z;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) {
	int N = result.size();

	Pixel bMinusA = b - a;
	vec3 step;
	step.x = float(bMinusA.x);
	step.y = float(bMinusA.y);
	step.z = bMinusA.zinv;
	step = step / myMax(N - 1.0f, 1.0f);

	vec3 current;
	current.x = float(a.x);
	current.y = float(a.y);
	current.z = a.zinv;

	for (int i = 0; i<N; ++i) {
		result[i].x = int(current.x);
		result[i].y = int(current.y);
		result[i].zinv = current.z;

		if (b.x < a.x) { // descending
			if (result[i].x < b.x)
				result[i].x = b.x;
		}
		else {
			if (result[i].x > b.x)
				result[i].x = b.x;
		}

		if (b.y < a.y) { // descending
			if (result[i].y < b.y)
				result[i].y = b.y;
		}
		else {
			if (result[i].y > b.y)
				result[i].y = b.y;
		}

		if (b.zinv < a.zinv) { // descending
			if (result[i].zinv < b.zinv)
				result[i].zinv = b.zinv;
		}
		else {
			if (result[i].zinv > b.zinv)
				result[i].zinv = b.zinv;
		}

		current = current + step;
	}
}

// For test 6 ==============================================================================================

void DrawWithIllumination() {
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	// Reset depthBuffer
	for (int i = 0; i < SCREEN_HEIGHT; i++)
		for (int j = 0; j < SCREEN_WIDTH; j++)
			depthBuffer[i][j] = 0.0f;

	for (size_t i = 0; i<triangles.size(); ++i) {
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		vertices[0].normal = triangles[i].normal;
		vertices[1].normal = triangles[i].normal;
		vertices[2].normal = triangles[i].normal;

		vertices[0].reflectance = triangles[i].color;
		vertices[1].reflectance = triangles[i].color;
		vertices[2].reflectance = triangles[i].color;

		DrawPolygonWithIllumination(vertices);
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void DrawPolygonWithIllumination(const vector<Vertex>& vertices) {
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<PixelI> projectedVertices(V);
	for (int i = 0; i < V; ++i) {
		VertexShader(vertices[i], projectedVertices[i]);
	}

	vector<PixelI> leftPixels;
	vector<PixelI> rightPixels;
	ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
	DrawColoredRows(leftPixels, rightPixels);
}

void ComputePolygonRows(const vector<PixelI>& vertexPixels, vector<PixelI>& leftPixels, vector<PixelI>& rightPixels) {
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();
	for (size_t i = 0; i < vertexPixels.size(); i++) {
		minY = myMin(minY, vertexPixels[i].y);
		maxY = myMax(maxY, vertexPixels[i].y);
	}
	int numOfRows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels = vector<PixelI>(numOfRows);
	rightPixels = vector<PixelI>(numOfRows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < numOfRows; i++) {
		leftPixels[i].x = numeric_limits<int>::max();
		leftPixels[i].y = minY + i;
		leftPixels[i].zinv = 0.0f;

		rightPixels[i].x = numeric_limits<int>::min();
		rightPixels[i].y = minY + i;
		leftPixels[i].zinv = 0.0f;
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (size_t i = 0; i < vertexPixels.size(); ++i) {
		int j = (i + 1) % vertexPixels.size(); // The next vertex

		PixelI delta = myAbs(vertexPixels[i] - vertexPixels[j]);
		int pixels = myMax(delta.x, delta.y) + 1;

		vector<PixelI> result(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], result);
		for (size_t k = 0; k < result.size(); k++) {
			int indexY = result[k].y - minY;
			if (leftPixels[indexY].x > result[k].x) {
				leftPixels[indexY] = result[k];
			}
			if (rightPixels[indexY].x < result[k].x) {
				rightPixels[indexY] = result[k];
			}
		}
	}
}

void DrawColoredRows(const vector<PixelI>& leftPixels, const vector<PixelI>& rightPixels) {
	for (size_t i = 0; i < leftPixels.size(); i++) {
		int rowLength = rightPixels[i].x - leftPixels[i].x + 1;

		vector<PixelI> result(rowLength);
		Interpolate(leftPixels[i], rightPixels[i], result);

		for (int j = 0; j < rowLength; j++) {
			PixelShader(result[j]);
		}
	}
}

void VertexShader(const Vertex& v, PixelI& p) {
	vec3 vprime = v.position - cameraPosition;
	vprime = vprime * getRotationMatrix();

	p.x = int(focalLength * vprime.x / vprime.z + SCREEN_WIDTH / 2.0f);
	p.y = int(focalLength * vprime.y / vprime.z + SCREEN_HEIGHT / 2.0f);
	p.zinv = 1.0f / vprime.z;

	// Compute illumination
	vec3 ncap = glm::normalize(v.normal);
	vec3 r = lightPos - v.position;
	vec3 rcap = glm::normalize(r);

	float intensityFactor = glm::dot(ncap, rcap);
	if (intensityFactor < 0.0f)
		intensityFactor = 0.0f;

	vec3 directLightPowerPerArea = (lightPower * intensityFactor) / (4 * PI * glm::length(r));

	p.illumination = v.reflectance * (directLightPowerPerArea + indirectLightPowerPerArea);
}

void PixelShader(const PixelI& p) {
	int x = p.x;
	int y = p.y;

	if (x < 0 || y < 0 || x >= SCREEN_WIDTH || y >= SCREEN_HEIGHT)
		return;

	if (p.zinv >= depthBuffer[x][y]) {
		depthBuffer[x][y] = p.zinv;
		PutPixelSDL(screen, x, y, p.illumination);
	}
}

void Interpolate(PixelI a, PixelI b, vector<PixelI>& result) {
	int N = result.size();

	PixelI bMinusA = b - a;
	vec3 step;
	step.x = float(bMinusA.x);
	step.y = float(bMinusA.y);
	step.z = bMinusA.zinv;
	step = step / myMax(N - 1.0f, 1.0f);

	vec3 stepI = bMinusA.illumination;
	stepI = stepI / myMax(N - 1.0f, 1.0f);

	vec3 current;
	current.x = float(a.x);
	current.y = float(a.y);
	current.z = a.zinv;

	vec3 currentI;
	currentI = a.illumination;

	for (int i = 0; i<N; ++i) {
		result[i].x = int(current.x);
		result[i].y = int(current.y);
		result[i].zinv = current.z;
		result[i].illumination = currentI;

		if (b.x < a.x) { // descending
			if (result[i].x < b.x)
				result[i].x = b.x;
		}
		else {
			if (result[i].x > b.x)
				result[i].x = b.x;
		}

		if (b.y < a.y) { // descending
			if (result[i].y < b.y)
				result[i].y = b.y;
		}
		else {
			if (result[i].y > b.y)
				result[i].y = b.y;
		}

		if (b.zinv < a.zinv) { // descending
			if (result[i].zinv < b.zinv)
				result[i].zinv = b.zinv;
		}
		else {
			if (result[i].zinv > b.zinv)
				result[i].zinv = b.zinv;
		}

		//TODO no bounds check for currentI (illumination) currently

		current = current + step;
		currentI = currentI + stepI;
	}
}

// For test 7/8 ==============================================================================================

void DrawWithPixelIllumination() {
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	// Reset depthBuffer
	for (int i = 0; i < SCREEN_HEIGHT; i++)
		for (int j = 0; j < SCREEN_WIDTH; j++)
			depthBuffer[i][j] = 0.0f;

	for (size_t i = 0; i<triangles.size(); ++i) {

		//// ** Start Backface-culling code **
		if (isBackface(triangles[i].v0, triangles[i].normal)) {
			continue;
		}
		//// ** End Backface-culling code **

		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;

		DrawPolygonWithPixelIllumination(vertices);
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void DrawPolygonWithPixelIllumination(const vector<Vertex>& vertices) {
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<PixelP> projectedVertices(V);
	for (int i = 0; i < V; ++i) {
		VertexShader(vertices[i], projectedVertices[i]);
	}

	vector<PixelP> leftPixels;
	vector<PixelP> rightPixels;
	ComputePolygonRows(projectedVertices, leftPixels, rightPixels);
	DrawColoredRows(leftPixels, rightPixels);
}

void ComputePolygonRows(const vector<PixelP>& vertexPixels, vector<PixelP>& leftPixels, vector<PixelP>& rightPixels) {
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int minY = numeric_limits<int>::max();
	int maxY = numeric_limits<int>::min();
	for (size_t i = 0; i < vertexPixels.size(); i++) {
		minY = myMin(minY, vertexPixels[i].y);
		maxY = myMax(maxY, vertexPixels[i].y);
	}
	int numOfRows = maxY - minY + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels = vector<PixelP>(numOfRows);
	rightPixels = vector<PixelP>(numOfRows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < numOfRows; i++) {
		leftPixels[i].x = numeric_limits<int>::max();
		leftPixels[i].y = minY + i;
		leftPixels[i].zinv = 0.0f;

		rightPixels[i].x = numeric_limits<int>::min();
		rightPixels[i].y = minY + i;
		leftPixels[i].zinv = 0.0f;
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (size_t i = 0; i < vertexPixels.size(); ++i) {
		int j = (i + 1) % vertexPixels.size(); // The next vertex

		PixelP delta = myAbs(vertexPixels[i] - vertexPixels[j]);
		int pixels = myMax(delta.x, delta.y) + 1;

		vector<PixelP> result(pixels);
		Interpolate(vertexPixels[i], vertexPixels[j], result);
		for (size_t k = 0; k < result.size(); k++) {
			int indexY = result[k].y - minY;
			if (leftPixels[indexY].x > result[k].x) {
				leftPixels[indexY] = result[k];
			}
			if (rightPixels[indexY].x < result[k].x) {
				rightPixels[indexY] = result[k];
			}
		}
	}
}

void DrawColoredRows(const vector<PixelP>& leftPixels, const vector<PixelP>& rightPixels) {
	for (size_t i = 0; i < leftPixels.size(); i++) {
		int rowLength = rightPixels[i].x - leftPixels[i].x + 1;

		vector<PixelP> result(rowLength);
		Interpolate(leftPixels[i], rightPixels[i], result);

		for (int j = 0; j < rowLength; j++) {
			PixelShader(result[j]);
		}
	}
}

void VertexShader(const Vertex& v, PixelP& p) {
	vec3 vprime = v.position - cameraPosition;
	vprime = vprime * getRotationMatrix();

	p.x = int(focalLength * vprime.x / vprime.z + SCREEN_WIDTH / 2.0f);
	p.y = int(focalLength * vprime.y / vprime.z + SCREEN_HEIGHT / 2.0f);
	p.zinv = 1.0f / vprime.z;
	p.pos3d = v.position;
}

void PixelShader(const PixelP& p) {
	int x = p.x;
	int y = p.y;

	if (x < 0 || y < 0 || x >= SCREEN_WIDTH || y >= SCREEN_HEIGHT)
		return;

	if (p.zinv >= depthBuffer[x][y]) {
		depthBuffer[x][y] = p.zinv;

		// ** Start Selective illumination **
		// do not compute directLight if backfacing
		//if (isBackface(p.pos3d, currentNormal)) {
		//	vec3 pixelColor = currentReflectance * (indirectLightPowerPerArea);
		//	PutPixelSDL(screen, x, y, pixelColor);
		//	return;
		//}
		// ** End Selective illumination **


		// Compute illumination
		vec3 ncap = glm::normalize(currentNormal);
		vec3 r = lightPos - p.pos3d;
		vec3 rcap = glm::normalize(r);

		float intensityFactor = glm::dot(ncap, rcap);
		if (intensityFactor < 0.0f)
			intensityFactor = 0.0f;

		vec3 directLightPowerPerArea = (lightPower * intensityFactor) / (4 * PI * glm::length(r));

		vec3 pixelColor = currentReflectance * (directLightPowerPerArea + indirectLightPowerPerArea);

		PutPixelSDL(screen, x, y, pixelColor);
	}
}

void Interpolate(PixelP a, PixelP b, vector<PixelP>& result) {
	int N = result.size();

	//a.pos3d = a.pos3d * a.zinv;
	//b.pos3d = b.pos3d * b.zinv;

	PixelP bMinusA = b - a;
	vec3 step;
	step.x = float(bMinusA.x);
	step.y = float(bMinusA.y);
	step.z = bMinusA.zinv;
	step = step / myMax(N - 1.0f, 1.0f);

	vec3 stepP = bMinusA.pos3d;
	stepP = stepP / myMax(N - 1.0f, 1.0f);

	vec3 current;
	current.x = float(a.x);
	current.y = float(a.y);
	current.z = a.zinv;

	vec3 currentP;
	currentP = a.pos3d;

	for (int i = 0; i<N; ++i) {
		result[i].x = int(current.x);
		result[i].y = int(current.y);
		result[i].zinv = current.z;
		result[i].pos3d = currentP;// / current.z;

		if (b.x < a.x) { // descending
			if (result[i].x < b.x)
				result[i].x = b.x;
		}
		else {
			if (result[i].x > b.x)
				result[i].x = b.x;
		}

		if (b.y < a.y) { // descending
			if (result[i].y < b.y)
				result[i].y = b.y;
		}
		else {
			if (result[i].y > b.y)
				result[i].y = b.y;
		}

		if (b.zinv < a.zinv) { // descending
			if (result[i].zinv < b.zinv)
				result[i].zinv = b.zinv;
		}
		else {
			if (result[i].zinv > b.zinv)
				result[i].zinv = b.zinv;
		}

		//TODO no bounds check for currentP (pos) currently

		current = current + step;
		currentP = currentP + stepP;
	}
}

// Supporting functions ==================================================================================================

Pixel operator-(Pixel const& lhs, Pixel const& rhs) {
	Pixel result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	result.zinv = lhs.zinv - rhs.zinv;
	return result;
}

PixelI operator-(PixelI const& lhs, PixelI const& rhs) {
	PixelI result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	result.zinv = lhs.zinv - rhs.zinv;
	result.illumination = lhs.illumination - rhs.illumination;
	return result;
}

PixelP operator-(PixelP const& lhs, PixelP const& rhs) {
	PixelP result;
	result.x = lhs.x - rhs.x;
	result.y = lhs.y - rhs.y;
	result.zinv = lhs.zinv - rhs.zinv;
	result.pos3d = lhs.pos3d - rhs.pos3d;
	return result;
}

template<class T> T myMax(T a, T b) {
	return (a > b) ? a : b;
}

template<class T> T myMin(T a, T b) {
	return (a < b) ? a : b;
}

PixelP myAbs(PixelP const& p) {
	PixelP result;
	result.x = myAbs(p.x);
	result.y = myAbs(p.y);
	result.zinv = myAbs(p.zinv);
	result.pos3d = glm::abs(p.pos3d);
	return result;
}

PixelI myAbs(PixelI const& p) {
	PixelI result;
	result.x = myAbs(p.x);
	result.y = myAbs(p.y);
	result.zinv = myAbs(p.zinv);
	result.illumination = glm::abs(p.illumination);
	return result;
}

Pixel myAbs(Pixel const& p) {
	Pixel result;
	result.x = myAbs(p.x);
	result.y = myAbs(p.y);
	result.zinv = myAbs(p.zinv);
	return result;
}

float myAbs(float const& v) {
	return (v < 0.0f) ? -v : v;
}

int myAbs(int const& v) {
	return (v < 0) ? -v : v;
}

// Back Face culling ==================================================================================================

// Check if a triangle is backwards facing in relation to the camera.
// Part of the Backface-culling procedure that is supposed to speed up rendering time.
// Assumptions: Triangle that is passed into this algorithm must not be a Bézier surface, 
//				ie, all points on the triangle must have the same normal.
bool isBackface(vec3 pos3d, vec3 triangleNormal) {
	vec3 cameraFacing = pos3d - cameraPosition;
	float dotProduct = glm::dot(cameraFacing, triangleNormal);

	if (dotProduct < 0) {
		return true;
	}
	else {
		return false;
	}
}

// Custom Loader ==================================================================================================

void LoadCustomModel(vector<Triangle>& triangles) {

	cout << "Loading custom model.." << endl;

	FILE * input;
	fopen_s(&input, "bun1.model.txt", "r");

	// Load points
	int pointsCount;
	fscanf_s(input, "%d", &pointsCount);

	vector<float> intensities;
	vector<vec3> points;
	points.reserve(pointsCount);


	float xf, yf, zf, intensity;
	for (int i = 0; i < pointsCount; i++) {
		fscanf_s(input, "%f %f %f %*f %f", &xf, &yf, &zf, &intensity);

		intensities.push_back(intensity);
		points.push_back(vec3(xf, -yf, -zf));
	}

	// Load triangles
	int trianglesCount;
	fscanf_s(input, "%d", &trianglesCount);

	triangles.clear();
	triangles.reserve(trianglesCount);

	vec3 color(1, 1, 1);

	int xi, yi, zi;
	for (int i = 0; i < trianglesCount; i++) {
		fscanf_s(input, "%*d %d %d %d", &xi, &yi, &zi);

		triangles.push_back(Triangle(points[xi], points[yi], points[zi], color));
	}

	cout << "Custom model successfully loaded!" << endl;
}

int LoadPlyFile(const char * filePath, vector<Triangle>& triangles, bool reverseX, bool reverseY, bool reverseZ) {
	FILE * input;
	fopen_s(&input, filePath, "rb");

	int mode = 0; // 1 = ASCII, 2 = binary_little_endian, 3 = binary_big_endian
	bool readVertices = false;
	int vertexCount = 0;
	bool readFaces = false;
	int faceCount = 0;

	int xyzCount = 0;

	vector<int> dataType; // 1 = signed integral, 2 = unsigned integral, 3 = floating-point
	vector<int> dataLength; // length of the data type, in bytes
	vector<int> dataIgnore; // 1 = ignore the data (unknown property), 0 = read it!

	vector<int> faceDataType; //same as above, but with 4 = the list.
	vector<int> faceDataLength;
	vector<int> faceDataIgnore;

	char str[100];
	char * token;

	cout << "Reading PLY file header" << endl;

	// === Read PLY identifier ===
	fgets(str, 99, input);
	if (strcmp(str, "ply\n") != 0) {
		return PLY_LOADER_INVALID_IDENTIFIER;
	}

	// === Read format ===
	token = readHeaderLine(input);
	if (token == NULL || strcmp(token, "format") != 0) {
		return PLY_LOADER_INVALID_FORMATTING;
	}
	token = strtok(NULL, " \t\r\n");
	if (token == NULL) {
		return PLY_LOADER_INVALID_FORMATTING;
	}
	if (strcmp(token, "ascii") == 0) {
		mode = 1;
	}
	else if (strcmp(token, "binary_little_endian") == 0) {
		mode = 2;
	}
	else if (strcmp(token, "binary_big_endian") == 0) {
		mode = 3;
	}
	else {
		return PLY_LOADER_UNKNOWN_FORMAT;
	}
	//ignore version

	// === Read each element ===
	token = readHeaderLine(input);
	while (true) {
		if (token == NULL) {
			return PLY_LOADER_INVALID_FORMATTING;
		}
		if (strcmp(token, "end_header") == 0) {
			break;
		}
		if (strcmp(token, "element") != 0) {
			return PLY_LOADER_INVALID_FORMATTING;
		}

		token = strtok(NULL, " \t\r\n");
		if (token == NULL) {
			return PLY_LOADER_INVALID_FORMATTING;
		}
		if (strcmp(token, "vertex") == 0) {
			readVertices = true;
			//get size
			token = strtok(NULL, " \t\r\n");
			if (token == NULL) {
				return PLY_LOADER_INVALID_FORMATTING;
			}
			sscanf(token, "%d", &vertexCount);

			//read vertex properties
			while (true) {
				token = readHeaderLine(input);
				if (token == NULL) {
					return PLY_LOADER_INVALID_FORMATTING;
				}
				if (strcmp(token, "property") == 0) {
					token = strtok(NULL, " \t\r\n");
					if (token == NULL) {
						return PLY_LOADER_INVALID_FORMATTING;
					}

					if (strcmp(token, "char") == 0) {
						dataType.push_back(1);
						dataLength.push_back(1);
					}
					else if (strcmp(token, "uchar") == 0) {
						dataType.push_back(2);
						dataLength.push_back(1);
					}
					else if (strcmp(token, "short") == 0) {
						dataType.push_back(1);
						dataLength.push_back(2);
					}
					else if (strcmp(token, "ushort") == 0) {
						dataType.push_back(2);
						dataLength.push_back(2);
					}
					else if (strcmp(token, "int") == 0) {
						dataType.push_back(1);
						dataLength.push_back(4);
					}
					else if (strcmp(token, "uint") == 0) {
						dataType.push_back(2);
						dataLength.push_back(4);
					}
					else if (strcmp(token, "float") == 0) {
						dataType.push_back(3);
						dataLength.push_back(4);
					}
					else if (strcmp(token, "double") == 0) {
						dataType.push_back(3);
						dataLength.push_back(8);
					}
					else {
						return PLY_LOADER_INVALID_FORMATTING;
					}

					token = strtok(NULL, " \t\r\n");
					if (token == NULL) {
						return PLY_LOADER_INVALID_FORMATTING;
					}

					if (strcmp(token, "x") == 0 || strcmp(token, "y") == 0 || strcmp(token, "z") == 0) {
						xyzCount++;
						dataIgnore.push_back(0);
					}
					else {
						dataIgnore.push_back(1);
					}
				}
				else {
					break;
				}
			}

			if (xyzCount != 3) {
				return PLY_LOADER_WRONG_NUMBER_OF_AXIS;
			}
		}
		else if (strcmp(token, "face") == 0) {
			readFaces = true;
			//get size
			token = strtok(NULL, " \t\r\n");
			if (token == NULL) {
				return PLY_LOADER_INVALID_FORMATTING;
			}
			sscanf(token, "%d", &faceCount);

			//read faces properties
			while (true) {
				token = readHeaderLine(input);
				if (token == NULL) {
					return PLY_LOADER_INVALID_FORMATTING;
				}
				if (strcmp(token, "property") == 0) {
					token = strtok(NULL, "\r\n");
					if (token == NULL) {
						return PLY_LOADER_INVALID_FORMATTING;
					}
					if (strcmp(token, "list uchar int vertex_indices") != 0) {
						return PLY_LOADER_UNSUPPORTED_FACE_PROPERTY;
					}
				}
				else {
					break;
				}
			}
		}
		else {
			//ignore the element, as it is unknown
			while (true) {
				token = readHeaderLine(input);
				if (token == NULL) {
					return PLY_LOADER_INVALID_FORMATTING;
				}
				if (strcmp(token, "property") == 0) {
					//continue;
				}
				else {
					break;
				}
			}
		}
	}

	// === Read the body ===
	if (mode == 1) { //ASCII
		// Construct the fscanf format string
		// In ASCII mode, everything is read as a double, despite its original data type, as scanf can hand the variance
		char readFormatString[1000] = "";
		for (int i = 0; i < dataType.size(); i++) {
			if (dataIgnore[i] == 1) {
				strcat(readFormatString, "%*lf ");
			}
			else {
				strcat(readFormatString, "%lf ");
			}
		}
		readFormatString[strlen(readFormatString) - 1] = '\0'; // remove last space

		// Load points
		vector<vec3> vertices;
		vertices.reserve(vertexCount);

		double xf, yf, zf;
		for (int i = 0; i < vertexCount; i++) {
			if (i % 100000 == 0) {
				cout << "Reading vertex " << i << " of " << vertexCount << endl;
			}

			fscanf_s(input, readFormatString, &xf, &yf, &zf);

			vertices.push_back(
				vec3(
					(reverseX) ? -xf : xf,
					(reverseY) ? -yf : yf,
					(reverseZ) ? -zf : zf
				)
			);
		}

		// Load triangles
		triangles.clear();
		triangles.reserve(faceCount);

		vec3 color(1, 1, 1);

		int ci, xi, yi, zi;
		for (int i = 0; i < faceCount; i++) {
			if (i % 100000 == 0) {
				cout << "Reading triangle " << i << " of " << faceCount << endl;
			}

			fscanf_s(input, "%d %d %d %d", &ci, &xi, &yi, &zi);

			if (ci != 3) {
				triangles.clear();
				return PLY_LOADER_UNSUPPORTED_POLYGON;
			}

			triangles.push_back(Triangle(vertices[xi], vertices[yi], vertices[zi], color));
		}
	}
	else { // 2 or 3, i.e. binary_little_endian or binary_big_endian
		// Construct the fscanf format string
		int dataPosition[] = { 0, 0, 0 }; //holds a reference to which index the x, y, z property is
		int dataPositionCounter = 0;

		char readFormatString[1000] = "";
		char tempStr[10];
		for (int i = 0; i < dataType.size(); i++) {
			if (dataIgnore[i] == 1) {
				strcat(readFormatString, "%*");
			}
			else {
				strcat(readFormatString, "%");
				dataPosition[dataPositionCounter++] = i;
			}

			sprintf(tempStr, "%d", dataLength[i]);

			strcat(readFormatString, tempStr);
			strcat(readFormatString, "c");
		}

		// Load points
		vector<vec3> vertices;
		vertices.reserve(vertexCount);

		unsigned char binaryData[3][8];
		float coord[3];

		for (int i = 0; i < vertexCount; i++) {
			if (i % 100000 == 0) {
				cout << "Reading vertex " << i << " of " << vertexCount << endl;
			}

			fscanf(input, readFormatString, binaryData[0], binaryData[1], binaryData[2]);

			for (int j = 0; j < 3; j++) {
				int k = dataPosition[j];
				if ((mode == 2 && systemIsBigEndian()) || (mode == 3 && !systemIsBigEndian())) {
					reverseEndianness(binaryData[j], dataLength[k]);
				}

				if (dataType[k] == 1 && dataLength[k] == 1) {
					coord[j] = float(*((signed char *)binaryData[j]));
				}
				else if(dataType[k] == 2 && dataLength[k] == 1) {
					coord[j] = float(*((unsigned char *)binaryData[j]));
				}
				else if (dataType[k] == 1 && dataLength[k] == 2) {
					coord[j] = float(*((signed short *)binaryData[j]));
				}
				else if (dataType[k] == 2 && dataLength[k] == 2) {
					coord[j] = float(*((unsigned short *)binaryData[j]));
				}
				else if (dataType[k] == 1 && dataLength[k] == 4) {
					coord[j] = float(*((signed int *)binaryData[j]));
				}
				else if (dataType[k] == 2 && dataLength[k] == 4) {
					coord[j] = float(*((unsigned int *)binaryData[j]));
				}
				else if (dataType[k] == 3 && dataLength[k] == 4) {
					coord[j] = *((float *)binaryData[j]);
				}
				else if (dataType[k] == 3 && dataLength[k] == 8) {
					coord[j] = float(*((double *)binaryData[j]));
				}
			}

			vertices.push_back(
				vec3(
					(reverseX) ? -coord[0] : coord[0],
					(reverseY) ? -coord[1] : coord[1],
					(reverseZ) ? -coord[2] : coord[2]
				)
			);
		}

		// Load triangles
		triangles.clear();
		triangles.reserve(faceCount);

		vec3 color(1, 1, 1);

		unsigned char ci;
		unsigned char binaryData2[8];
		int indexes[3];
		for (int i = 0; i < faceCount; i++) {
			if (i % 100000 == 0) {
				cout << "Reading triangle " << i << " of " << faceCount << endl;
			}

			fscanf(input, "%1c%4c%4c%4c", &ci, binaryData[0], binaryData[1], binaryData[2]);

			if (ci != unsigned char(3)) {
				triangles.clear();
				return PLY_LOADER_UNSUPPORTED_POLYGON;
			}

			for (int j = 0; j < 3; j++) {
				if ((mode == 2 && systemIsBigEndian()) || (mode == 3 && !systemIsBigEndian())) {
					reverseEndianness(binaryData[j], 4);
				}
				
				indexes[j] = *((int *)binaryData[j]);
			}

			triangles.push_back(Triangle(vertices[indexes[0]], vertices[indexes[1]], vertices[indexes[2]], color));
		}
	}

	cout << "PLY file loaded." << endl;

	return PLY_LOADER_OK;
}

char * readHeaderLine(FILE * input) {
	char str[1000];
	char * token;

	while (true) {
		fgets(str, 999, input);

		token = strtok(str, " \t\r\n");

		if (token == NULL) {
			return NULL;
		}
		if (strcmp(token, "comment") != 0) {
			return token;
		}
	}
}

int systemIsBigEndian() {
	union {
		uint32_t i;
		char c[4];
	} bint = { 0x01020304 };

	return bint.c[0] == 1;
}

void reverseEndianness(unsigned char * binaryData, int size) {
	char temp;

	for (int i = 0; i < size / 2; i++) {
		temp = binaryData[i];
		binaryData[i] = binaryData[size - i - 1];
		binaryData[size - i - 1] = temp;
	}
}

void printBounds(vector<Triangle>& triangles) {
	float mins[3] = { numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max() };
	float maxs[3] = { numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min() };
	for (int i = 0; i < triangles.size(); i++) {
		mins[0] = myMin(mins[0], triangles[i].v0.x);
		mins[0] = myMin(mins[0], triangles[i].v1.x);
		mins[0] = myMin(mins[0], triangles[i].v2.x);
		maxs[0] = myMax(mins[0], triangles[i].v0.x);
		maxs[0] = myMax(mins[0], triangles[i].v1.x);
		maxs[0] = myMax(mins[0], triangles[i].v2.x);

		mins[1] = myMin(mins[1], triangles[i].v0.y);
		mins[1] = myMin(mins[1], triangles[i].v1.y);
		mins[1] = myMin(mins[1], triangles[i].v2.y);
		maxs[1] = myMax(mins[1], triangles[i].v0.y);
		maxs[1] = myMax(mins[1], triangles[i].v1.y);
		maxs[1] = myMax(mins[1], triangles[i].v2.y);

		mins[2] = myMin(mins[2], triangles[i].v0.z);
		mins[2] = myMin(mins[2], triangles[i].v1.z);
		mins[2] = myMin(mins[2], triangles[i].v2.z);
		maxs[2] = myMax(mins[2], triangles[i].v0.z);
		maxs[2] = myMax(mins[2], triangles[i].v1.z);
		maxs[2] = myMax(mins[2], triangles[i].v2.z);
	}

	cout << mins[0] << " <= x <= " << maxs[0] << endl;
	cout << mins[1] << " <= y <= " << maxs[1] << endl;
	cout << mins[2] << " <= z <= " << maxs[2] << endl;
}