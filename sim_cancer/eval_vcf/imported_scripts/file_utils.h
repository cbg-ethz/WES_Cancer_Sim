#ifndef __FILE_UTILS_H__
#define __FILE_UTILS_H__

#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>

int my_mkdir(char* dirname)
{
	char command[10000];
	sprintf(command, "mkdir -p %s", dirname);
	int ret = system(command); 

	if (ret!=0)
	{
		fprintf(stderr, "[%s] Error: command '%s' failed: %i\n", __func__, command, ret); 
	}
	return ret; 
}

#define IS_OTHER 0
#define IS_DIR 1
#define IS_REG_FILE 3
#define IS_SLINK 2

int file_stats(char* filename)
{
	struct stat sb;

	if (stat(filename, &sb) == -1)
	{
		return -1;
	}
	
	switch (sb.st_mode & S_IFMT) 
	{
		case S_IFBLK:	return IS_OTHER;	break;//block device
		case S_IFCHR:	return IS_OTHER;	break;//character device
		case S_IFDIR:	return IS_DIR;	break;//directory
		case S_IFIFO:	return IS_OTHER;	break;//FIFO/pipe
		case S_IFLNK:	return IS_SLINK;	break;//symlink
		case S_IFREG:	return IS_REG_FILE;	break;//regular file
		case S_IFSOCK:	return IS_OTHER;	break;//socket
		default:		return IS_OTHER;	break;
	}
}

bool fexist(char* filename)
{
	return file_stats(filename)>IS_DIR; // symlink or regular file
}
#endif
