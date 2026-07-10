#pragma once


#include <gl/glsl.h>
#include <gl/texture.h>
#include <gl/vertexbuffer.h>
#include <gl/indexbuffer.h>

#include <vector> 


struct Node
{
	Node() {}
    virtual ~Node() {};
};

struct Child;
struct Parent 
{
    std::vector<Child*> children;
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
	// Geometry borrows its buffers. Their owner must outlive the geometry node.
	VertexBuffer*	vbo;
	IndexBuffer*	ibo;

    Geometry() : vbo(NULL), ibo(NULL) {}
    virtual ~Geometry() {};
};

// TODO typedef int RenderTarget;
struct NodeState;
struct SceneState;
struct Renderable
{
    vec4 AABB;
    virtual void Render(NodeState& nodeState, SceneState& sceneState ) = 0;
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
struct Updatable : Buffer
{
    virtual void Update() = 0;
    virtual ~Updatable() {};
};



