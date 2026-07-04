#pragma once

#include <map>
#include <vector>
#include <assert.h>

namespace gd
{
    typedef unsigned long long Entity;
/*    typedef void (*EntityDestroyerFn)(Entity);

    inline std::vector<EntityDestroyerFn>& GetEntityDestroyers()
    {
        static std::vector<EntityDestroyerFn> destroyers;
        return destroyers;
    }
*/
    Entity CreateEntity()
    {
        static Entity nextID = 0;
        return nextID++;
    }

    // Replacement for typeid(T) or std::type_index 
    template<typename T>
    size_t GetTypeID()
    {
        static size_t id = 0;
        return id++;
    }






    // 
    template <typename T>
    std::map<Entity, std::vector<T>>& GetComponentsMap();

    template <typename T>
    void RemoveComponents(Entity entity)
    {
        GetComponentsMap<T>().erase(entity);
    }
/*
    template <typename T>
    struct RegisterComponentDestroyer
    {
        RegisterComponentDestroyer()
        {
            GetEntityDestroyers().push_back(&RemoveComponents<T>);
        }
    };
*/
    template <typename T>
    std::map<Entity, std::vector<T>>& GetComponentsMap()
    {
        //static RegisterComponentDestroyer<T> registerDestroyer;
        static std::map<Entity, std::vector<T>> componentsMap;
        return componentsMap;
    }
/*
    void DestroyEntity(Entity entity)
    {
        std::vector<EntityDestroyerFn>& destroyers = GetEntityDestroyers();
        for (size_t i = 0; i < destroyers.size(); ++i)
        {
            destroyers[i](entity);
        }
    }
*/
    template <typename T>
    std::vector<T>& GetComponents(Entity node)
    {
        std::map<Entity, std::vector<T>>& componentsMap = GetComponentsMap<T>();
        std::vector<T>& components = componentsMap[node];
        return components;
    }

    template <typename T>
    T& AddComponent(Entity node)
    {
        std::vector<T>& components = GetComponents<T>(node);
        T component;
        components.push_back(component);
        return components.back();
    }
/*
    // TODO Parent
    // TODO Children
    struct _Hierarchy
    {
	    Entity parent;
	    std::vector<Entity> children;
    };

    Entity AddChild(Entity node)
    {
        std::vector<_Hierarchy>& components = GetComponents<_Hierarchy>(node);
        _Hierarchy& nodeHierarchy = (components.size() == 0) ? AddComponent<_Hierarchy>(node) : components[0];

        Entity child = CreateEntity();
        _Hierarchy& childHierarchy = AddComponent<_Hierarchy>(child);

        nodeHierarchy.children.push_back(child);
        childHierarchy.parent = node;

        return child;
    }

    // TODO GetChildren();
*/

};




/*
void main_ECS()
{
    typedef float _vec3[3];

    struct _Transform
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


    Entity node = CreateEntity();

    //
    // Assert AddComponent creates a component for the node and can be later on retrieved.
    //

    std::vector<_Transform>& components0 = GetComponents<_Transform>(node);
    assert(components0.size() == 0);

    _Transform& t0 = AddComponent<_Transform>(node);
    assert(components0.size() == 1);

    t0.position[0] = 123;

    std::vector<_Transform>& components1 = GetComponents<_Transform>(node);
    assert(components1.size() == 1);
    assert(&components0 == &components1);

    _Transform& t1 = components1.back();

    assert(t1.position[0] = 123);
    assert(&t0 == &t1);

    //
    // Scene example with hierarchy
    //

    Entity scene = CreateEntity();
    _Transform& camera = AddComponent<_Transform>(scene);

    assert(GetComponents<_Transform>(scene).size() == 1);

    Entity heightmap = AddChild(scene);
    _Transform& heightmapTransform = AddComponent<_Transform>(heightmap);

	std::vector<_Hierarchy>& sceneHierarchyPool = GetComponents<_Hierarchy>(scene);
    assert(sceneHierarchyPool.size() == 1);
    assert(sceneHierarchyPool[0].children.size() == 1);
    assert(sceneHierarchyPool[0].children[0] == heightmap);

    std::vector<_Hierarchy>& heightmapHierarchyPool = GetComponents<_Hierarchy>(heightmap);
    assert(heightmapHierarchyPool.size() == 1);
    assert(heightmapHierarchyPool[0].parent = scene);
}
*/