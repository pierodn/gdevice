#pragma once

#include "type/glsl.h"
#include "type/array.h"

// NOTE
//   Allocation/deallocation to video memory
//   
//   Allocation to video memory is performed by Renderer on first attempt of using the resource.
//   
//   Deallocation is performed by the resource itself on destruction
//     by calling the abstract "deallocator" function, applied to the given video memory id
//
//   If invalidated then it means the copy in video memory requires to be updated.


struct Cacheable
{
	uint id;
	void (*deallocator)(uint&);

	bool videomem_invalidated;
	bool hostmem_invalidated;

    Cacheable() : id(0), deallocator(NULL), videomem_invalidated(true), hostmem_invalidated(true)
	{
	}

    ~Cacheable()
    {
		deallocate();
    }

	void deallocate()
	{
		if( deallocator && id ) deallocator(id);
	}
};
