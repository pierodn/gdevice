#pragma once

//#include <windows.h>
#include "os/platform.h"
/*
#ifdef _OPENMP

#include <omp.h>

struct Mutex
{
	omp_lock_t olock;

	Mutex() { omp_init_lock(&olock); }
	~Mutex() { omp_destroy_lock(&olock); }
	void lock() { omp_set_lock(&olock); }
	void unlock() { omp_unset_lock(&olock); }
};

#endif
*/

#if 1

class Mutex
{
	CRITICAL_SECTION critical_section;
public:

	Mutex()
	{
		InitializeCriticalSection( &critical_section );
	}

	~Mutex()
	{
		DeleteCriticalSection( &critical_section );
	}

	void lock()
	{
		EnterCriticalSection( &critical_section );
	}

	bool tryLock()
	{
		return 0 != TryEnterCriticalSection( &critical_section );
	}

	void unlock()
	{
		LeaveCriticalSection( &critical_section );
	}
};

#elif 0

class Mutex
{
	bool critical_section;

public:

	Mutex()
	{
		critical_section = false;
	}

	~Mutex()
	{
	}

	void lock()
	{
		while( critical_section );
		critical_section = true;
	}

	void unlock()
	{
		while( !critical_section );
		critical_section = false;
	}
};

#elif 0

class Mutex
{
	CRITICAL_SECTION critical_section;
/*
	void print_data()
	{
		printf("LockCount       = %i\n", critical_section.LockCount );
		printf("LockSemaphore   = %x\n", critical_section.LockSemaphore );
		printf("OwningThread    = %x\n", critical_section.OwningThread );
		printf("RecursionCount  = %i\n", critical_section.RecursionCount );
		printf("SpinCount       = %i\n", critical_section.SpinCount );
		printf("\n");
	}
*/
public:

	Mutex()
	{
		//printf("---- create mutex ----\n");
		InitializeCriticalSection( &critical_section );
		//print_data();
	}

	~Mutex()
	{
		//printf("---- destroy mutex ----\n");
		DeleteCriticalSection( &critical_section );
		//print_data();
	}

	void lock()
	{
		//printf("---- locking ----\n");
		//print_data();
		//if( critical_section.LockCount<0 )
			EnterCriticalSection( &critical_section );
		//print_data();
	}

	void unlock()
	{
		//printf("---- unlocking ----\n");
		//print_data();
		//if( critical_section.LockCount>=0 )
			LeaveCriticalSection( &critical_section );
		//print_data();
	}
};

#elif 0

class Mutex
{
	HANDLE m_hMutex;

public:
	enum eMutexState
	{
		eStopped,
		eMutexAquired,
		eTimedOut,
	};

	Mutex()
	{
		m_hMutex = CreateMutex(NULL, FALSE, NULL);

		if( m_hMutex == INVALID_HANDLE_VALUE )
		{
			printf("Mutex(): m_hMutex == INVALID_HANDLE_VALUE !!\n" );
		}
	}

	~Mutex()
	{
		if (m_hMutex != INVALID_HANDLE_VALUE)
			CloseHandle(m_hMutex);
	}
 
	bool IsValid() const             { return m_hMutex != INVALID_HANDLE_VALUE; }
/*
	eMutexState lock(HANDLE hStopEvent, DWORD dwTimeout = INFINITE)
	{
		ASSERT(IsValid());

		HANDLE hArray[2] = { hStopEvent, m_hMutex };

		switch (WaitForMultipleObjects(2, hArray, FALSE, dwTimeout))
		{
			case WAIT_OBJECT_0:
				return eStopped;

			case WAIT_OBJECT_0 + 1:
				return eMutexAquired;
		
			default:
				return eTimedOut;
		}
	}
	*/
	void lock( HANDLE hStopEvent )
	{
		ASSERT(IsValid());

		HANDLE hArray[2] = { hStopEvent, m_hMutex };

		switch (WaitForMultipleObjects(2, hArray, FALSE, INFINITE))
		{
			case WAIT_OBJECT_0:
				return eStopped;

			case WAIT_OBJECT_0 + 1:
				return eMutexAquired;
		
			default:
				return eTimedOut;
		}
	}

	void unlock() const        { ::ReleaseMutex(m_hMutex); }

};

#endif
