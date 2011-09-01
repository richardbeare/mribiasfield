#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "Debug.h"

#ifdef WIN32
	#undef  NDEBUG
	#include <assert.h>
#endif

#define tesss


int TAB_LEVEL = 0;
const int TAB_SPACE_COUNT = 2;


void DebugHelper::DUMP(const char* formatStr, ...)
{
	va_list		argList;
	char		outBuffer[1024];

	int spaceCount = TAB_LEVEL * TAB_SPACE_COUNT;
	for(int i = 0; i < spaceCount; i++) 
	{
		outBuffer[i] = ' ';
	}

	va_start(argList, formatStr);
	vsprintf(outBuffer + spaceCount, formatStr, argList);
	va_end(argList);

	printf("%s", outBuffer);
	fflush(stdout);
}

void DebugHelper::DUMP_AND_EXIT(const char* formatStr, ...)
{
	va_list		argList;
	char		outBuffer[1024];

	va_start(argList, formatStr);
	vsprintf(outBuffer, formatStr, argList);
	va_end(argList);

	printf("%s\n", outBuffer);
	fflush(stdout);

#ifdef WIN32
	assert(0);
#endif
	//printf("Press enter to continue.\n", outBuffer);
	//fflush(stdout);
	//getchar();	

	exit(EXIT_FAILURE);
}


