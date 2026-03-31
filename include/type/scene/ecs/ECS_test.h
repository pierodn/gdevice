#include <stdio.h>


#if 1
// #define _DEBUG

//#define TSL_NO_EXCEPTIONS
#define _HAS_EXCEPTIONS 0
//#define _DEFINE_DEPRECATED_HASH_CLASSES 0
//#include <tr1/unordered_map>
//#include <hash_map>
#include <map>
#include <vector>
#include <assert.h>

// TODO gd::Entity and gd::Node (which is ECS Entity with Inheritance).
// TODO rename Component as DataPart, Attribute, Property ?
// TODO Apply logic to components of a node (update function).


//namespace gd
//{
    typedef unsigned long Node2;

    Node2 AddNode()
    {
        static Node2 nextID = 0;
        return ++nextID;
    }

    template <typename T>
    std::map<Node2, std::vector<T>>& GetComponentsMap()
    {
        static std::map<Node2, std::vector<T>> componentsMap;
        return componentsMap;
    }

    template <typename T>
    std::vector<T>& GetComponents(Node2 node)
    {
        std::map<Node2, std::vector<T>>& componentsMap = GetComponentsMap<T>();
        std::vector<T>& components = componentsMap[node];
        return components;
    }

    template <typename T>
    T& AddComponent(Node2 node)
    {
        std::vector<T>& components = GetComponents<T>(node);
        T component;
        components.push_back(component);
        return components.back();
    }


    struct Hierarchy
    {
	    Node2 parent;
	    std::vector<Node2> children;
    };

    Node2 AddChild(Node2 node)
    {
        std::vector<Hierarchy>& components = GetComponents<Hierarchy>(node);
        Hierarchy& nodeHierarchy = (components.size() == 0) ? AddComponent<Hierarchy>(node) : components[0];

        Node2 child = AddNode();
        Hierarchy& childHierarchy = AddComponent<Hierarchy>(child);

        nodeHierarchy.children.push_back(child);
        childHierarchy.parent = node;

        return child;
    }

    // TODO GetChildren();

//};




typedef float _vec3[3];

struct Transform
{
    _vec3 position;
    _vec3 rotation;
    _vec3 scale;
};


struct SomeDataPart
{
    int count;
    char* base;
    long count2;
};

void main_ECS()
{
    Node2 node = AddNode();

    //
    // Assert AddComponent creates a component for the node and can be later on retrieved.
    //

    std::vector<Transform>& components0 = GetComponents<Transform>(node);
    assert(components0.size() == 0);

    Transform& t0 = AddComponent<Transform>(node);
    assert(components0.size() == 1);

    t0.position[0] = 123;

    std::vector<Transform>& components1 = GetComponents<Transform>(node);
    assert(components1.size() == 1);
    assert(&components0 == &components1);

    Transform& t1 = components1.back();

    assert(t1.position[0] = 123);
    assert(&t0 == &t1);

    //
    // Scene example with hierarchy
    //

    Node2 scene = AddNode();
    Transform& camera = AddComponent<Transform>(scene);

    assert(GetComponents<Transform>(scene).size() == 1);

    Node2 heightmap = AddChild(scene);
    Transform& heightmapTransform = AddComponent<Transform>(heightmap);

	std::vector<Hierarchy>& sceneHierarchyPool = GetComponents<Hierarchy>(scene);
    assert(sceneHierarchyPool.size() == 1);
    assert(sceneHierarchyPool[0].children.size() == 1);
    assert(sceneHierarchyPool[0].children[0] == heightmap);

    std::vector<Hierarchy>& heightmapHierarchyPool = GetComponents<Hierarchy>(heightmap);
    assert(heightmapHierarchyPool.size() == 1);
    assert(heightmapHierarchyPool[0].parent = scene);
}

#endif

void main()
{
    main_ECS();
    printf("hello spank\n");
}