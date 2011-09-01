#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <stdio.h>
#include <string>

class DebugHelper
{
public:
	static void DUMP(const char *formatStr, ...);
	static void DUMP_AND_EXIT(const char *formatStr, ...);
};

extern int TAB_LEVEL;
extern const int TAB_SPACE_COUNT;

#define TOUCH(x)        ;{x;};

#define PRINT_STR(x)    printf(#x ": %s\n",   (char *)       (x));
#define PRINT_STRING(x) printf(#x ": %s\n",   (char *)       (x).c_str());
#define PRINT_INT(x)    printf(#x ": %i\n",   (int)          (x));
#define PRINT_TS(x)     printf(#x ": %lld\n", (long long int)(x));
#define PRINT_UINT(x)   printf(#x ": %u\n",   (unsigned int) (x));
#define PRINT_FLOAT(x)  printf(#x ": %f\n",   (float)        (x));
#define PRINT_SCI(x)    printf(#x ": %e\n",   (double)       (x));
#define PRINT_DOUBLE(x) printf(#x ": %f\n",   (double)       (x));
#define PRINT_OBJ(x)    {                       \
                           printf(#x ":\n");    \
						   x.Dump();            \
                        }

#define PRINT_INTS(x)   {                                                              \
							printf(#x ": [ ");                                         \
							for(int i = 0; i < (x).size(); i++) printf("%i ", (x)[i]); \
							printf("]\n");                                             \
                        }

#define REPORT_LOC()      printf("[%s, %i]: \n", __FILE__, __LINE__);
#define DEBUG_DUMP	      REPORT_LOC(); DebugHelper::DUMP
#define REPORT_FAILURE    REPORT_LOC(); DebugHelper::DUMP_AND_EXIT
#define REPORT            DebugHelper::DUMP

#define INCR_FUNC_LEVEL()  {                                                \
						 	 REPORT("\nExecuting %s()\n\n", __FUNCTION__);  \
							 TAB_LEVEL++;                                   \
                           }

#define DCR_FUNC_LEVEL()   {                          \
                              REPORT("\nDone.\n");    \
						  	  TAB_LEVEL--;            \
							  INSIST(TAB_LEVEL >= 0); \
                           }

#define INCR_TAB_LEVEL()   {                                            \
							 TAB_LEVEL++;                               \
                           }

#define DCR_TAB_LEVEL()    {                          \
						  	  TAB_LEVEL--;            \
							  INSIST(TAB_LEVEL >= 0); \
                           }


#define ENSURE(expression)   if((expression) == false)                            \
							 {                                                    \
								 REPORT_FAILURE("ENSURE: %s\n", #expression);     \
							 }

#define INSIST(expression)   if((expression) == false)                            \
							 {                                                    \
							     REPORT_FAILURE("INSIST: %s\n", #expression);     \
							 }

#define CHECK_RANGE(x, min, max) {                            \
								    ENSURE((x) >= (min));     \
									ENSURE((x) <= (max));     \
							   	 }

#define INSIST_RANGE(x, min, max) {                           \
								    INSIST((x) >= (min));     \
									INSIST((x) <= (max));     \
							   	  }

#endif


