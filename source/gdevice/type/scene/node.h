#pragma once


#include "type/glsl.h"
#include "type/array.h"

#include "__temp/texture.h"
#include "gpu/program.h"
#include "__temp/vertexbuffer.h"
#include "gpu/indexbuffer.h"

#include "gpu/renderer.h"

struct Node
{
	Node() {}
    virtual ~Node() {};
};

struct Child;
struct Parent 
{
	Array<Child*> children; // TODO use Child* instead of Node*
    virtual ~Parent() {};
};

struct Child 
{
    Parent* parent;

    Child() : parent(NULL) {};
    virtual ~Child() {};
};

struct Transform 
{
    mat3 transform;

    Transform()
	{
		transform.position	= vec3(0);
		transform.rotation	= vec3(0);
		transform.scale		= vec3(1);
    }
    virtual ~Transform() {};
};

struct Geometry 
{
	// TODO Material*		material; // TerrainMaterial by default
    // TODO Abstract vertex attribute => so to allow other configurations.
    // TODO VAO => stores the state of vertex attributes for performance.
	VertexBuffer*	vbo;
	IndexBuffer*	ibo;

    Geometry() : vbo(NULL), ibo(NULL) {}
    virtual ~Geometry() 
    {
        // TODO delete vbo, ibo
    };
};

typedef int RenderTarget;
class Renderer;
struct Renderable
{
    vec4 AABB;
    virtual void Render(Renderer& renderer, RenderTarget& target) = 0;
    virtual ~Renderable() {};
};

struct Impostor 
{
    vec2			impostorDirection;
	float			impostorDistance;
	Texture<rgba>	impostorTexture; // TODO Updatable
    virtual ~Impostor() {};
};

class Renderer;
struct Updatable : Cacheable
{
    virtual void Update(Renderer& renderer) = 0;
    virtual ~Updatable() {};
};



