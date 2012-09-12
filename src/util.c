#include<string.h>
#include<R.h>
#include"Mutils.h"

/* Lapack condition number approximation: currently only supports _1 or _Inf norm : */
char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error("argument type[1]='%s' must be a character string of string length 1",
	      typstr);
    typup = toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error("argument type[1]='%s' must be one of '1','O', or 'I'",
	      typstr);
    return typup; /* 'O' or 'I' */
}
