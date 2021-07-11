#pragma once

#include <windows.h>



class Thread
{
	HANDLE handle;

public:

	Thread()
	{
		handle = NULL;
	}

	~Thread()
	{
		terminate();
	}

	void run( void (*function)(void*), void far* data, int core=0, int priority=0 )
	{
		if( !function ) return;

		handle = CreateThread( NULL, 0, (PTHREAD_START_ROUTINE) function, (LPVOID)data, 0, NULL );
	}

	inline static void sleep( unsigned long milliseconds )
	{
		Sleep( milliseconds );
	}

	void suspend()
	{
		SuspendThread( handle );
	}

	void resume()
	{
		ResumeThread( handle );
	}

	void terminate()
	{
		if( handle && TerminateThread(handle, 0) )
		{
			CloseHandle( handle );
		}
		handle = NULL;
	}

	bool isAlive()
	{
		return handle>0/* && WaitForSingleObject( handle, 0 )*/;
	}
/*
	void wait( unsigned long milliseconds = INFINITE )
	{
		WaitForSingleObject( handle, milliseconds );
	}
*/	



};