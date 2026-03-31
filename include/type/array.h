#pragma once

template<class T>
struct Array
{
	T* array;
	int length;
	int allocated;
	int granularity;

	Array(int allocated = 0, int granularity = 1) : 
        array(NULL), 
        length(0), 
        allocated(allocated), 
        granularity(granularity) 
	{
        if( allocated ) 
        {
            allocate(allocated);
        }
	}

	~Array() 
	{ 
		delete[] array; 
	}
/*
	// copy constructor
	// TODO test this against memory leaks
	void operator=( const Array& copy )
	{
		array = NULL;
		granularity = copy.granularity;
		length = copy.length;
		allocate( length );
		for(int i=0; i<length; i++) array[i] = copy.array[i];
	}
*/
	const T& operator[]( int index ) const 
	{
		return array[index]; 
	}

	T& operator[]( int index ) 
	{
		return array[index]; 
	}

	Array& push(const T& value)
	{
        if(length + 1 > allocated) {
            allocate(allocated+granularity);
        }
		array[length++] = value;
		return *this;
	}

	const T pop()
	{
		return array[--length];
	}

	Array& push()
	{
        if(length + 1 > allocated) {
            allocate(allocated + granularity);
        }
		length++;
		return *this;
	}

	T& tail()
	{
		return array[length-1]; 
	}

	Array& allocate(int allocated, int granularity = 0)
	{
		if(this->allocated != allocated) {
			array = (T*) realloc(array, allocated*sizeof(T));
			this->allocated = allocated;
		}

		if(granularity>0) {
			this->granularity = granularity;
		}

		return *this;
	}

	Array& fit()
	{
		return allocate(length);
	}

	Array& reset()
	{
		delete[] array;
		array = NULL;
		length = 0;
		allocated = 0; 
		//granularity = 1; //
		return *this;
	}

	Array& clean()
	{
		memset(array, 0, allocated*sizeof(T));
		return *this;
	}

	int size()
	{
		return length;
	}
};
